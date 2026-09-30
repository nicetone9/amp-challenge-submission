"""Sequence-only VQ codec, causal code prior, and official DiMA components."""
from types import SimpleNamespace
import torch
from torch import nn
import torch.nn.functional as F
from .vendor.models.score_estimator import ScoreEstimator
from .vendor.diffusion.schedulers import Tanh
from .vendor.diffusion.dynamic import DynamicSDE
from .vendor.diffusion.solvers import DDIMSolver

AA = "ACDEFGHIKLMNPQRSTVWY"

class Residual(nn.Module):
    def __init__(self):
        super().__init__()
        self.c1 = nn.Conv1d(128, 128, 3, padding=1)
        self.c2 = nn.Conv1d(128, 128, 3, padding=1)
    def forward(self, x, mask):
        return (x + self.c2(F.gelu(self.c1(x)) * mask)) * mask

class Codec(nn.Module):
    def __init__(self):
        super().__init__()
        self.input = nn.Linear(1280, 128)
        self.blocks = nn.ModuleList([Residual(), Residual()])
        self.output = nn.Linear(128, 64)
        self.decoder = nn.Sequential(nn.Linear(64,128),nn.GELU(),nn.Linear(128,20))
        self.register_buffer("book", torch.zeros(512,64))
        self.register_buffer("ema_count", torch.ones(512))
        self.register_buffer("ema_sum", torch.zeros(512,64))
        self.register_buffer("window", torch.zeros(512,dtype=torch.long))
    def encode(self, x, mask):
        h = self.input(x).transpose(1,2) * mask[:,None]
        for block in self.blocks:
            h = block(h, mask[:,None])
        return self.output(h.transpose(1,2)) * mask[...,None]
    def forward(self, x, mask):
        z = self.encode(x,mask)
        with torch.no_grad():
            flat = z[mask]
            distance = flat.square().sum(-1,keepdim=True) + self.book.square().sum(-1) - 2*flat@self.book.T
            ids = torch.zeros(mask.shape,dtype=torch.long,device=x.device)
            ids[mask] = distance.argmin(-1)
            q = self.book[ids] * mask[...,None]
        logits = self.decoder(z + (q-z).detach())
        commitment = ((z-q.detach()).square().mean(-1)*mask).sum(-1)/mask.sum(-1)
        return logits, commitment, ids, z
    @torch.no_grad()
    def ema_update(self, ids, latents, step):
        counts = torch.bincount(ids,minlength=512).float()
        sums = torch.zeros_like(self.ema_sum).index_add_(0,ids,latents)
        self.ema_count.mul_(.99).add_(counts,alpha=.01)
        self.ema_sum.mul_(.99).add_(sums,alpha=.01)
        self.book.copy_(self.ema_sum/self.ema_count[:,None].clamp_min(1e-6))
        self.window.add_(counts.long())
        reset = 0
        if step % 1000 == 0:
            dead = self.window.eq(0)
            reset = int(dead.sum())
            if reset:
                replacements = latents[torch.randint(len(latents),(reset,),device=latents.device)]
                self.book[dead] = replacements
                self.ema_sum[dead] = replacements
                self.ema_count[dead] = 1
            self.window.zero_()
        return reset

class Prior(nn.Module):
    def __init__(self):
        super().__init__()
        self.embedding = nn.Embedding(514,256,padding_idx=513)
        self.position = nn.Embedding(50,256)
        self.length = nn.Embedding(51,256)
        layer = nn.TransformerEncoderLayer(256,8,1024,dropout=.1,batch_first=True,norm_first=True)
        self.transformer = nn.TransformerEncoder(layer,6,norm=nn.LayerNorm(256),enable_nested_tensor=False)
        self.output = nn.Linear(256,512)
    def forward(self, tokens, lengths):
        n = tokens.shape[1]
        mask = torch.arange(n,device=tokens.device)[None] < lengths[:,None]
        x = self.embedding(tokens)+self.position(torch.arange(n,device=tokens.device))[None]+self.length(lengths)[:,None]
        causal = torch.ones(n,n,dtype=torch.bool,device=tokens.device).triu(1)
        return self.output(self.transformer(x,mask=causal,src_key_padding_mask=~mask))

class Decoder(nn.Module):
    def __init__(self, config):
        super().__init__()
        from transformers import EsmConfig
        from transformers.models.esm.modeling_esm import EsmLMHead
        self.head = EsmLMHead(EsmConfig(**config))
        # Replaced by the checkpoint's exact tokenizer-derived indices.
        self.register_buffer("aa_ids", torch.zeros(20, dtype=torch.long))
    def forward(self,x):
        return self.head(x).index_select(-1,self.aa_ids)

def denoiser():
    return ScoreEstimator(SimpleNamespace(
        embedding_size=1280,hidden_size=320,num_hidden_layers=12,
        num_attention_heads=16,attention_head_size=20,intermediate_size=1280,
        max_position_embeddings=50,use_self_cond=True,qk_norm=True,
        attention_probs_dropout_prob=0.,hidden_dropout_prob=0.,
        layer_norm_eps=1e-12,add_cross_attention=False))

def dynamic():
    return DynamicSDE(Tanh(10),T=1)

def sequence_ce(logits,target,mask):
    loss = F.cross_entropy(logits.transpose(1,2),target,ignore_index=-100,reduction="none")
    return (loss*mask).sum(-1)/mask.sum(-1)
