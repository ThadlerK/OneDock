# app/utils.py
import yaml
import os
import streamlit as st

CONFIG_PATH = "config/config.yaml"

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