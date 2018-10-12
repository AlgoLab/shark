# Read-filter
Tool for mapping-free classification of reads, using Bloom filter.

## Dependencies
Read-filter requires the following libraries and tools:

- [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)
- [KMC v3](https://github.com/refresh-bio/KMC/tree/v3.0.0)
- [htslib v1.8](https://github.com/samtools/htslib/tree/1.8)

This repository comes with them as submodules.

## Download and Installation
In order to install the repository run the followings steps.
First clone the repository and move into it.

```shell
git clone --recursive https://github.com/tmunimib/read-filter.git
cd read-filter
```

If you have KMC3, sdsl-lite, and htslib already installed you can skip the following commands.
Please Note: KMC seems to not compile using g++-7 under Ubuntu; please switch to g++-6 in this step to avoid compilation errors.

```shell
chmod +x Library_installer.sh
./Library_installer.sh
```

Now you can compile read-filter from the root of you local copy of the repository.
```shell
cd <path-to-read-filter-local-repo>
make
```    
