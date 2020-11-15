import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cozmic", # Replace with your own username
    version="0.0.1",
    author=["Lluís Revilla Sancho","Ferran Muiños"],
    author_email=["lluis.revilla@gmail.com","ferran.muinos@gmail.com"],
    description="Correlated mutations and distance correlations to predict aminoacid interactions.",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/llrs/PYT-SBI/",
    packages=["cozmic"],
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.4",
    install_requires=[
        "numpy",
        "matplotlib",
        "Biopython"
#        "logging",
#        "argparse",
#        "copy",
#        "distutils",
#        "ftplib",
#        "math",
#        "os",
#        "urllib",
    ]
)
