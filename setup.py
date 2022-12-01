"""
GUD: Genomic Unification Database
"""

import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="GUD",
    version="22.12.1",
    author="Oriol Fornes, Tamar Av-Shalom",
    author_email="oriol.fornes@gmail.com, avshalom.tamar0@gmail.com",
    description="Genomic Unification Database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wassermanlab/GUD",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
