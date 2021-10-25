"""Setup for package asdeep."""
from setuptools import setup

des = "A deep-learning tool to interpret variant function by allelic imbalance."
setup(
    name="asdeep",
    version="0.1.3",
    author="Zhenhua Zhang",
    author_email="zhenhua.zhang217@gmail.com",
    maintainer="Zhenhua Zhang",
    maintainer_email="zhenhua.zhang217@gmail.com",
    description=des,
    packages=["asdeep"],
    package_data={
        "asdeep": ["test.fa.gz", "test.vcf.gz", "test.gtf.gz", "test.bed.gz"]
    },
    platform="any"
)
