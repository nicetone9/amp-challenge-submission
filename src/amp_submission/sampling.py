"""Sampling copied from the evaluated 2026-09-29 pipeline, without training imports."""
import json
import os
import random
import numpy as np
import torch
from .models import AA, Codec, Prior, Decoder, denoiser, dynamic
from .vendor.diffusion.solvers import DDPMSolver

OBJECTIVES = ("broad", "gram_positive", "gram_negative", "mdr", "selectivity", "joint")

def seed_all(seed):
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    torch.backends.mha.set_fastpath_enabled(False)
    torch.use_deterministic_algorithms(True)

def load_models(root, arch, device):
    def restore(model, name):
        state = torch.load(root / name, map_location="cpu", weights_only=True)
        model.load_state_dict(state, strict=True)
        return model.to(device).eval().requires_grad_(False)
    policies = {o: restore(Prior() if arch == "vq" else denoiser(),
                           f"{arch}/{o}.pt") for o in OBJECTIVES}
    if arch == "vq":
        bundle = {"codec": restore(Codec(), "vq/codec.pt")}
        reference = None
    else:
        reference = restore(denoiser(), "dima/reference.pt")
        config = json.loads((root / "esm-config.json").read_text())
        norm = json.loads((root / "normalization.json").read_text())
        bundle = {"decoder": restore(Decoder(config), "dima/decoder.pt"),
                  "mean": torch.tensor(norm["mean"], device=device),
                  "std": torch.tensor(norm["std"], device=device)}
    return policies, reference, bundle

@torch.no_grad()
def sample(arch, policy, reference, bundle, lengths, device):
    lengths = torch.as_tensor(lengths, device=device, dtype=torch.long)
    width = int(lengths.max())
    if arch == "vq":
        inputs = torch.full((len(lengths), width), 513, device=device, dtype=torch.long)
        inputs[:, 0] = 512
        codes = torch.zeros_like(inputs)
        for position in range(width):
            logits = policy(inputs[:, :position + 1], lengths)[:, -1]
            codes[:, position] = torch.multinomial(torch.softmax(logits, -1), 1).squeeze(-1)
            if position + 1 < width:
                inputs[:, position + 1] = codes[:, position]
        aa = bundle["codec"].decoder(bundle["codec"].book[codes]).argmax(-1)
    else:
        mask = torch.arange(width, device=device)[None] < lengths[:, None]
        xt = torch.randn(len(lengths), width, 1280, device=device) * mask[..., None]
        conditioning = torch.zeros_like(xt)
        times = torch.linspace(.999, .001, 51, device=device)
        solver = DDPMSolver(dynamic(), None)
        for i in range(50):
            t, nt = times[i].expand(len(lengths)), times[i + 1].expand(len(lengths))
            conditioning = reference(xt, t, mask.float(), conditioning) * mask[..., None]
            x0 = policy(xt, t, mask.float(), conditioning) * mask[..., None]
            mu, sigma = solver.q_x_t_reverse(xt, x0, t, nt)
            xt = (mu + sigma.clamp_min(1e-8) * torch.randn_like(xt)) * mask[..., None]
        aa = bundle["decoder"](xt * bundle["std"] + bundle["mean"]).argmax(-1)
    return ["".join(AA[v] for v in row[:int(n)].cpu().tolist()) for row, n in zip(aa, lengths)]
