import setuptools

setuptools.setup(
    name="scdemon",
    version="0.1.0",
    url="https://github.edu/KellisLab/scdemon",
    author="Carles Boix",
    author_email="cboix@mit.edu",
    description="Computes and plots gene modules for single-cell data",
    long_description=open("README.md").read(),
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "scanpy",
        "anndata",
        "fbpca",
        "gprofiler-official",
        "python-igraph",
        "leidenalg",
        "adjustText",
        "umap-learn",
        "numba",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
    package_data={"": ["*.txt", "*.tsv"]},
)
