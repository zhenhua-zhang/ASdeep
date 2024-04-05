#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import torch
import pyro
import pyro.distributions as dist
import pyro.distributions.constraints as const

from pyro.infer import SVI
from pyro.infer import Trace_ELBO
from pyro.infer.mcmc import NUTS
from pyro.infer.mcmc.api import MCMC
from pyro.optim import Adam


out_dir = Path("~/Documents/git/ASdeep").expanduser()
n_samples = 1000
n_test = 1000

pyro.clear_param_store()
pyro.set_rng_seed(3141592)

def model(n, k):
    alpha = pyro.param("alpha", torch.tensor(1.0), constraint=const.positive)
    beta = pyro.param("beta", torch.tensor(1.0), constraint=const.positive)

    with pyro.plate("data", len(n)):
        pyro.sample("obs", dist.BetaBinomial(alpha, beta, n), obs=k) # Obervations

def guide(n, k):
    alpha = pyro.param("alpha", torch.tensor(1.0))
    beta = pyro.param("beta", torch.tensor(1.0))

n_vec = torch.randint(100, 300, (n_samples, ), dtype=torch.int32)
k_vec = torch.tensor([int(x * 0.95) for x in n_vec], dtype=torch.int32)

# Why with render_model the inferring results matches expectation?
pyro.render_model(model, model_args=(n_vec, k_vec), filename=out_dir/"temp/pyro_model.png", render_distributions=True)

stragegy = "svi"
if stragegy == "mcmc":
    nuts_kernel = NUTS(model, jit_compile=False, step_size=1e-5)
    MCMC(nuts_kernel, num_samples=1000, warmup_steps=1000, num_chains=1).run(n_vec, k_vec)
elif stragegy == "svi":
    adam_params = {"lr": 0.001, "betas": (0.95, 0.999)}
    optimizer = Adam(adam_params)
    svi = SVI(model, guide, optimizer, loss=Trace_ELBO())
    for _ in range(5000):
        svi.step(n=n_vec, k=k_vec)

alpha_hat = float(pyro.param("alpha").item())
beta_hat = float(pyro.param("beta").item())

bb = dist.BetaBinomial(alpha_hat, beta_hat, n_test)
print(f"alpha_hat is {alpha_hat:.4} and beta_hat is {beta_hat:.4}. BetaBinomial mean is {bb.mean:.4}")

x_val = np.linspace(1, n_test, n_test)
pdf_val = bb.log_prob(torch.tensor(x_val)).exp().numpy()

fig, axe = plt.subplots(1, 1)
_ = axe.plot(x_val, pdf_val)
fig.savefig(out_dir / "temp/pyro_test.png")
