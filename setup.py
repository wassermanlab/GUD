import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="GUD",
    version="0.0.1",
    author="Oriol Fornes",
    author_email="oriol.fornes@gmail.com",
    description="The Genomic Universal Database package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oriolfornes/GUD",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)