import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="joanne",  # Replace with your own username
    version="0.0.1",
    author="Geet George",
    author_email="geet.george@mpimet.mpg.de",
    description="EUREC4A Dropsonde Dataset",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Geet-George/JOANNE",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
