Docker file for the analysis of Mouse Gastruloids.

**Warning**: The analysis has not been completed yet. Some parts of the code ar nonsense so be careful to run it for now.

**Warning**: The bash files have been writen for Linux. Maybe in Mac or PC, the lines are different. The scripts are pretty simple so they should be straighforward to adapt.

# Use without Docker

In the file `requirements.txt` there is the defined basic packages used and its versions for the analysis. The analsys is performed using Python3.8.8.

Inside home, you will find all the analysis scripts.

The necessary data for the analysis can be obtained running

```
./download_data
```

# Use with Docker

All the analysis has been done in a docker image for reproducibility purposes. Check it in [dockerhub](https://hub.docker.com/r/dsblab/single_cell_analysis) or [github](https://github.com/dsb-lab/Docker-single-cell-analysis/tree/v0.2).

For running the docker and go over the analysis steps in a jupyter lab session, just run,

```
./docker/run.sh
```

Once run, the docker will be working and executing in port 8888. Change the channel in the `run.sh` script if you want to run it in some other channel. In the terminal you will find a token code that is generated for security reasons:

![](assets/token.png)

Copy that number. For accessing the session, open your favorite folder brwoser and search `localhost:8888`. Directly from the terminal it will be something like,

```
firefox localhost::8888
```

it will open a jupyterlab session that will ask for a password or token.  

![](assets/jupyterlab.png)

Copy the token you obtained before and you will be set up! 

Anything in the `home` folder can be seen by the docker and you and will be saved after the docker has finished.