from distutils.core import setup

setup(
    name="CANTHERM",
    version="0.1.0",
    author="James E T Smith",
    author_email="james.smith9113@gmail.com",
    packages=["cantherm"],
    license="LICENSE",
    description="Thermochemistry post processing for quantum chemistry packages",
    long_description=open("README.md").read(),
    python_requires=">3.5",
    install_requires=["scipy", "numpy", "matplotlib", "cclib>=1.6.3", "pytest-cov"],
)
