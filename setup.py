import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="joanne",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
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
