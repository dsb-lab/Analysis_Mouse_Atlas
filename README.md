Docker file for the analysis of Mouse Gastruloids.

**Warning**: The analysis has not been completed yet. Some parts of the code ar nonsense so be careful to run it for now.

**Warning**: The bash files have been writen for Linux. Maybe in Mac or PC, the lines are different. The scripts are pretty simple so they should be straighforward to adapt.

# Repository structure

The repository is structured as follows:

 - **Analysis**: Contains the scripts of the analysis for each dataset individually. So far, the following analysis are present:
    2.1. Pijuan: Preprocessing of the Pijuan dataset from the raw count matrix until the annotation of the different clusters of the analysis.
 - **assets**: Auxiliar images for the this README file.
 - **download_data.sh**: Bash file that automatizes the process of downloading the raw data and storing it in the Data folder.
 - **docker_run.sh**: Bash file that launches the docker in a jupyter lab environment for working with the analysis.

# Setting up the environment

## Download the data automatically

In order to prepare for a fresh start of the analysis it is required that the present folder contains a folder called Data. This folder can be link folder that redirects to any place where you want to store the data as long as you can access this place with your account privileges. For downloading automatically the data, just execute the commant:

```
./download_data.sh
```

This can take a while.

## Download the data by hand (if the previous did not work)

If the script doens't work for you. You can download the data by hand at the following places:

Pijuan dataset:
> https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz

When extracting the data from the `atlas_data.tar.gz` file, you will have to have, at least, the following data:

> Data/Pijuan/raw/barcodes.tsv
> Data/Pijuan/raw/genes.tsv
> Data/Pijuan/raw/raw_counts.mtx
> Data/Pijuan/raw/meta.tab

# Execute the jupyter lab environment

In order to go over the hole analysis with reproducibility, we developed a docker environment for this situation with all the tools already installed.

For running the environment, you will require to install [Docker](https://www.docker.com/) in your computer. Once done, just execute in a terminal inside the path of the repository:

```
./docker_run.sh
```

Once run, the docker will be working and executing in port the port specified by the variable `CHANNEL` in the docker_run.sh file. Change the channel in the  script if you want to run it in some other channel. In the terminal you will find a token code that is generated for security reasons:

![](assets/token.png)

Copy that number. For accessing the session, open your favorite folder brwoser and search `localhost:8888`. Directly from the terminal it will be something like,

```
firefox localhost::8888
```
where the 8888 if the channel you chose. It will open a jupyterlab session that will ask for a password or token.  

![](assets/jupyterlab.png)

Copy the token you obtained before and you will be prompted to the analysis environment! 
