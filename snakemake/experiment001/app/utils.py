# app/utils.py
import yaml
import os
import streamlit as st
import shutil

CONFIG_PATH = "config/config.yaml"
DATA_DIRS = [
    "data/inputs",
    "data/inputs/library_split",
    "data/interim", 
    "data/results",
]

def save_config(new_data):
    """
    Updates the config.yaml file with new keys/values 
    without overwriting existing ones.
    """
    os.makedirs("config", exist_ok=True)
    
    # Load existing config if it exists
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        config = {}

    # Update with new data
    config.update(new_data)

    # Write back to file
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(config, f)
    
    return config

def reset_project():
    """
    Wipes all data and resets configuration.
    """
    # 1. Delete the entire data folder
    if os.path.exists("data"):
        try:
            shutil.rmtree("data")
        except Exception as e:
            return f"Error deleting data: {e}"

    # 2. Re-create empty structure
    for d in DATA_DIRS:
        os.makedirs(d, exist_ok=True)

    # 3. Reset Config File to empty
    with open(CONFIG_PATH, "w") as f:
        yaml.dump({}, f)

    return "Success"