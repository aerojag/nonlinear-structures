# Standard library imports
import os, sys, json, subprocess, re
from pathlib import Path
from getpass import getuser
from string import Template

# Third party imports
import jinja2
import yaml
import dynamic_yaml

# -----------------------
# Global variables
# -----------------------
project_folder = os.getcwd()
base_dir = str(Path(os.path.dirname(os.path.realpath(__file__))).parents[0])  # Get path to the root folder of this project
res = f'{base_dir}/resources' # Path to the resources folder
docs = f'{base_dir}/docs'
user_id = getuser()