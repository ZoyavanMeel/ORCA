from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="orca",
    version="1.0",
    description="Tool for predicting the origin of replication on circular bacterial chromosomes",
    package_dir={"": "src"},
    packages=find_packages(where="src", exclude=["test", "Machine_learning", "data"]),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ZoyavanMeel/ORCA",
    author="Zoya van Meel",
    author_email="zoyavanmeel@gmail.com",
    license="GPL-3.0",
    classifiers=[
        "License :: GPL-3.0 License",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy >= 1.22.3",
        "scikit-learn >= 1.2.0",
        "biopython >= 1.79",
        "scipy >= 1.8.0",
        "matplotlib >= 3.5.2"
    ],
    extras_require={
        "dev": ["twine>=4.0.2"],
    },
    python_requires=">=3.10",
)
