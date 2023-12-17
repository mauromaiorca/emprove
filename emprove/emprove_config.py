#!/usr/bin/python3
import os
import json
from pathlib import Path
import argparse
import shutil


def setup_config(relion_reconstruct_path=None):
    user_config_dir = Path.home() / '.emprove'
    user_config_file = user_config_dir / 'config.json'
    user_config_dir.mkdir(exist_ok=True)
    config = {}
    if relion_reconstruct_path is not None:
        config["relion_reconstruct_path"] = relion_reconstruct_path
        config["relion_reconstruct_path_mpi"] = relion_reconstruct_path+"_mpi"  
    else:
        relion_reconstruct_path = shutil.which("relion_reconstruct")
        config["relion_reconstruct_path"] = relion_reconstruct_path
        config["relion_reconstruct_path_mpi"] = relion_reconstruct_path+"_mpi"          
        print ("relion_reconstruct_path=",relion_reconstruct_path)
    with open(user_config_file, 'w') as f:
        json.dump(config, f, indent=4)


def load_config():
    user_config_file = Path.home() / '.emprove' / 'config.json'
    if user_config_file.exists():
        with open(user_config_file, 'r') as f:
            config = json.load(f)
        return config, user_config_file
    else:
        return None, None


def main():
    parser = argparse.ArgumentParser(description="""
    This is a configuration script for the emprove library. 
    emprove is a library for assessing EM reconstruction from particle scores. 

    This script allows users to specify and change the path to the relion_reconstruct executable 
    which is an important part of the emprove setup. If the path is not provided, the script will 
    attempt to locate the executable automatically. 

    It also allows the user to overwrite existing configurations and see the current configuration. 
    The script stores the configuration in a JSON file in the user's home directory.
    """)
    parser.add_argument('--relion_reconstruct_path', help='Path to relion reconstruct')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing configuration')
    args = parser.parse_args()
    print ("args.relion_reconstruct_path=",args.relion_reconstruct_path)
    

    config, user_config_file = load_config()
    if config is not None:
        print(f"Configuration file location: {user_config_file}")
        print(f"Current Configuration:\n{json.dumps(config, indent=4)}")
        if args.overwrite and args.relion_reconstruct_path:
            print("Overwriting configuration.")
            setup_config(args.relion_reconstruct_path)
            print("New Configuration:\n{json.dumps(load_config()[0], indent=4)}")
    else:
        print("No configuration file found, creating file.")
        setup_config(args.relion_reconstruct_path)

if __name__ == "__main__":
    main()

