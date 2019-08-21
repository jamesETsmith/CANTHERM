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
    install_requires=["scipy", "numpy", "matplotlib", "cclib", "pytest-cov"],
)
