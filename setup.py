from setuptools import setup, find_packages

setup(
    name="q2-eco-processes",
    version="0.1.0",
    packages=find_packages(),
    author="LaBiOmics / UMC",
    author_email="bioinformatics@labiomics.org",
    description="Official QIIME 2 plugin for microbial ecological assembly processes quantification.",
    long_description=open("README.md").read() if open("README.md") else "",
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/LaBiOmics/q2-eco-processes",
    entry_points={
        'qiime2.plugins': ['q2-eco-processes=q2_eco_processes.plugin_setup:plugin'],
        'console_scripts': ['q2-eco-processes-check=q2_eco_processes._check:check_environment']
    },
    include_package_data=True,
    zip_safe=False
)
