
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# setuptools.setup(
#     name="fits_reader",
#     version="0.1.0",
#     author="RH",
#     author_email="kuma@kwansei.ac.jp",
#     description="fits read module'",
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url="https://github.com/RH-REP/CIBER2/fits_reader",
#     packages=setuptools.find_packages(),
#     classifiers=[
#         "Programming Language :: Python :: 3",
#         "License ::  Approved :: MIT License",
#         "Operating System :: OS Independent",
#     ],
#     # entry_points = {
#     #     'console_scripts': ['sample_command = sample_command.sample_command:main']
#     # },
#     python_requires='>=3.7',
# )

setuptools.setup(
    name="ciber2",
    version="0.1.0",
    author="RH",
    author_email="kuma@kwansei.ac.jp",
    description="spec read module'",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RH-REP/CIBER2",
#     packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License ::  Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=['spec_reader', 'fits_reader'],
    # entry_points = {
    #     'console_scripts': ['sample_command = sample_command.sample_command:main']
    # },
    python_requires='>=3.7',
)
