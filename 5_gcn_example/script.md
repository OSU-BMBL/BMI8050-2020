# Login OnDemand:

https://ondemand.osc.edu/pun/sys/dashboard

# Start a Jupyter session with the following settings:

Cluster: owens

Project: PAS1791

Hours: 4

Node type: any GPU

CUDA version: 10.1 **(important!!!)**

Cores: 4

Jupyter lab version: 2.1 (Default)

# Login OSC OnDemand and open a terminal session, the following command will update the course Git repository

```
cd ~/BMI8050-2020

git stash && git pull && git checkout stash -- .
```

# Install networkx

```
ml python/3.6
source activate bmi8050
pip install --user networkx
```