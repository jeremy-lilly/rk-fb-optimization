# Optimization of Runge-Kutta time-stepping schemes with forward-backward averaging

*Author: Jeremy Lilly*

This repo contains Jupyter notebooks written in Python for obtaining "optimal" time-stepping schemes for the shallow-water equations (SWEs).


## First-time setup

It is suggested to run these notebooks from within a virtual python environment (venv) rather than messing with the primary Python installation on your machine. To create a new venv:
```
python3 -m venv <path-to-new-venv>
```
For example, `python3 -m $HOME/my-venv` will create a new virtual python environment called `my-venv` in your home directory.

Next, activate the new venv and install the required packages with `pip`:
```
source <path-to-new-venv>/bin/activate
pip install -r requirements.txt
```

You can then leave the venv with:
```
deactivate
```

## Usage

To run these notebooks after the above first-time setup, simply activate the venv:
```
source <path-to-new-venv>/bin/activate
```
Open the Jupyter notebook interface:
```
jupyter notebook
```
Then, use the GUI to navigate to and run the notebooks.

To shutdown the notebook server, simply hit CTRL-C in the terminal where the server was launched and follow the prompts. To leave the venv, run:
```
deactivate
```

