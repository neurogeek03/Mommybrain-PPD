"""Quick test: check GPU visibility and key imports."""
import torch
print(f"PyTorch: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU count: {torch.cuda.device_count()}")
    print(f"GPU name: {torch.cuda.get_device_name(0)}")
    print(f"CUDA version: {torch.version.cuda}")
    # Quick tensor op on GPU
    x = torch.randn(1000, 1000, device="cuda")
    y = x @ x.T
    print(f"GPU matmul test: OK ({y.shape})")
else:
    print("No GPU detected!")

import scvi
print(f"scvi-tools: {scvi.__version__}")

import scanpy as sc
print(f"scanpy: {sc.__version__}")

try:
    import scib
    print(f"scib: {scib.__version__}")
except ImportError:
    print("scib: NOT INSTALLED")

print("\nAll imports OK.")
