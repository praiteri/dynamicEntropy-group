#!/usr/bin/env python3
import os
import sys
import glob
import yaml
import shutil
import logging
import subprocess
from colorama import Fore, Style
import random
import itertools

class ColoredFormatter(logging.Formatter):
    """
    Custom formatter that adds colors to log messages based on level with fixed width.
    """
    COLORS = {
        "DEBUG": Fore.CYAN + Style.BRIGHT,
        "INFO": Fore.WHITE,
        "WARNING": Fore.YELLOW + Style.BRIGHT,
        "ERROR": Fore.RED + Style.BRIGHT,
        "CRITICAL": Fore.MAGENTA + Style.BRIGHT,
    }
    LEVEL_NAME_WIDTH = 8

    def format(self, record):
        orig_levelname = record.levelname
        orig_msg = record.msg
        color_code = self.COLORS.get(orig_levelname, "")
        padding = " " * (self.LEVEL_NAME_WIDTH - len(orig_levelname))
        record.levelname = f"{color_code}{orig_levelname}{padding}{Style.RESET_ALL}"
        record.msg = f"{color_code}{orig_msg}{Style.RESET_ALL}"
        result = super().format(record)
        record.levelname = orig_levelname
        record.msg = orig_msg
        return result

class QuickScript:
    """
    Main class to handle script execution, folder management, and logging.
    """
    def __init__(self):
        self.dirs = ["./"]
        self.dryrun = False
        self.if_folder_exists = "enter"
        self.skip_if_exist = None
        self._batch_mode = False

        formatter = ColoredFormatter(fmt="%(levelname)8s - %(message)s")
        handler = logging.StreamHandler(stream=sys.stdout)
        handler.setFormatter(formatter)
        self.logger = logging.getLogger("script")
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.WARNING)

    def _handle_folders(self, folder):
        """
        Ensure the folder exists or handle it according to the specified behavior.
        """
        if not isinstance(folder, str):
            raise TypeError(f"Invalid type for 'dirs': {type(folder)}")
        if folder in ["./", "."]:
            return
        if os.path.isdir(folder):
            if self.if_folder_exists == "replace":
                self.logger.warning(f"Overwriting folder {folder}")
                shutil.rmtree(folder)
                os.makedirs(folder)
            elif self.if_folder_exists == "create":
                self.logger.error(f"Folder {folder} already exists")
                sys.exit(1)
            elif self.if_folder_exists == "skip":
                self.logger.warning(f"Folder {folder} already exists, skipping")
            elif self.if_folder_exists == "enter":
                pass
            else:
                self.logger.error(f"Unknown folder exists behavior: {self.if_folder_exists}")
                sys.exit(1)
        else:
            if self.if_folder_exists == "enter":
                self.logger.error(f"Folder {folder} doesn't exist!")
                sys.exit(1)
            else:
                os.makedirs(folder)
                self.logger.debug(f"Creating folder {folder}")

    def _check_already_processed(self):
        """
        Check if files in skip_if_exist have already been processed.
        """
        if self.skip_if_exist is None:
            return False
        already_processed = all(os.path.isfile(fp) for fp in self.skip_if_exist)
        if already_processed:
            self.logger.info("Folder already processed")
        return already_processed

    def _execute_commands(self, commands):
        """
        Execute a list of shell commands and return the results.
        """
        if isinstance(commands, str):
            commands = [commands]
        elif not isinstance(commands, list):
            raise ValueError("Commands should be a list of strings or a single string")
        if self.dryrun:
            for command in commands:
                self.logger.info(f"Running: {command}")
            return ["Dry run"] * len(commands)
        if self._check_already_processed():
            return ["Already processed"]
        results = []
        for command in commands:
            self.logger.info(f"Running: {command}")
            if command.startswith("cd "):
                os.chdir(command[3:])
                res = ""
            else:
                res = subprocess.run(
                    command, shell=True, text=True, capture_output=True, check=False
                )
                if res.returncode != 0:
                    self.logger.error("Command failed")
                    for line in res.stderr.strip().split("\n"):
                        self.logger.error(f"{line}")
                    sys.exit(1)
                else:
                    for line in res.stdout.strip().split("\n"):
                        self.logger.debug(f"{line}")
                    self.logger.debug("-" * 80)
            results.append(res)
        return results

    def _alive_bar(self, total=None):
        """
        Fallback implementation when alive_progress is not available.
        """
        class DummyBar:
            def __init__(self, total=None):
                self.total = total
            def __enter__(self):
                return self.__call__
            def __exit__(self, *args):
                pass
            def __call__(self, *args, **kwargs):
                pass
        if total is None:
            return DummyBar(total)
        try:
            from alive_progress import alive_bar, config_handler
            config_handler.set_global(force_tty=True, receipt=False, monitor=False)
            return alive_bar(total)
        except ImportError:
            return DummyBar(total)

    def dry_run(self, dryrun=True):
        self.dryrun = dryrun

    def set_noise_level(self, level):
        """
        Set the logging level.
        """
        levels = {"INFO": logging.INFO, "DEBUG": logging.DEBUG, "WARNING": logging.WARNING}
        if level.upper() in levels:
            self.logger.setLevel(levels[level.upper()])
        else:
            raise Exception(f"Unknown logging level: {level}")

    def folder_handler(self, behaviour):
        """
        Set behavior for existing folders.
        """
        if behaviour not in ["replace", "create", "skip", "enter"]:
            raise ValueError(f"Unknown behaviour for existing folders: {behaviour}")
        self.if_folder_exists = behaviour

    def find_matching_entry(self, dictionary, target_list):
        """
        Find a matching entry in a dictionary based on a target list.
        """
        for key, value in dictionary.items():
            try:
                clean_key = key.strip('()')
                parts = clean_key.split(',')
                parsed_items = []
                for part in parts:
                    part = part.strip()
                    if part.startswith('"') and part.endswith('"'):
                        parsed_items.append(part.strip('"'))
                    else:
                        try:
                            parsed_items.append(int(part))
                        except ValueError:
                            parsed_items.append(part)
                if len(parsed_items) != len(target_list):
                    self.logger.error(
                        f"\nList length mismatch:\nParsed items: {len(parsed_items)} elements {parsed_items}\nTarget list: {len(target_list)} elements {target_list}"
                    )
                    sys.exit(1)
                if all(parsed_items[i] == '*' or parsed_items[i] == target_list[i] for i in range(len(target_list))):
                    return value
            except Exception as e:
                self.logger.warning(e)
                continue
        self.logger.error(f"Cannot find variable combination that matches the entry {target_list}")
        sys.exit(1)

    def run_script(self, list_of_commands, dirs=None, variables=None):
        """
        Execute shell commands across multiple directories and collect their results.
        """
        original_dir = os.getcwd()
        result = {}
        if dirs is None:
            dirs = ["./"]
        if variables is None:
            variables = [None]
        ntmp = None
        if not self._batch_mode:
            ntmp = len(self.dirs) * len(variables)
        with self._alive_bar(ntmp) as bar:
            idir = -1
            for vv in variables:
                idir += 1
                for dd in dirs:
                    if vv is not None:
                        for k, v in vv.items():
                            dd = dd.replace(f"{{{{{k}}}}}", str(v))
                    os.chdir(original_dir)
                    self._handle_folders(dd)
                    if not os.path.isdir(dd):
                        self.logger.error(f"Directory {dd} does not exist!")
                    os.chdir(dd)
                    self.logger.info(f"Working in folder: {dd}")
                    for cmd in list_of_commands:
                        if isinstance(cmd, str):
                            lcmd = cmd
                        elif isinstance(cmd, list):
                            lcmd = cmd[idir]
                        elif isinstance(cmd, dict):
                            keys_list = list(vv.values())
                            lcmd = self.find_matching_entry(cmd, keys_list)
                        else:
                            self.logger.error(f"Unknown command type: {type(cmd)}")
                            sys.exit(1)
                        if vv is not None:
                            for k, v in vv.items():
                                lcmd = lcmd.replace(f"{{{{{k}}}}}", str(v))
                        lcmd = lcmd.replace("{{RANDOM}}", str(random.randint(0, 10000)))
                        result[dd] = self._execute_commands(lcmd)
                    bar()
                    self.logger.info("-" * 80)
                    os.chdir(original_dir)
        return result

def read_yaml_safe(file_path):
    """
    Safely read a YAML file and return its contents as a Python object.
    """
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            data = yaml.safe_load(file)
        return data
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def create_dict_combinations(input_dict):
    """
    Create all possible combinations of dictionary values.
    Converts numeric strings to integers automatically.
    """
    keys = list(input_dict.keys())
    values = list(input_dict.values())
    combinations = itertools.product(*values)
    result = []
    for combo in combinations:
        new_dict = {}
        for i, key in enumerate(keys):
            value = combo[i]
            try:
                new_dict[key] = int(value)
            except ValueError:
                new_dict[key] = value
        result.append(new_dict)
    return result

def main():
    if len(sys.argv) == 1:
        print("Usage: python run_script.py <yaml_file>")
        sys.exit(1)
    file_path = sys.argv[1]
    input_cmd = read_yaml_safe(file_path)
    qs = QuickScript()
    noise = input_cmd.get("noise", "INFO")
    qs.set_noise_level(noise)
    dry_run = input_cmd.get("dry_run", False)
    qs.dry_run(dry_run)
    qs._batch_mode = input_cmd.get("batch_mode", False)
    variables = input_cmd.get("variables", None)
    if variables is not None and not isinstance(variables, dict):
        raise TypeError(f"Invalid type for 'variables': {type(variables)}")
    list_of_variables = create_dict_combinations(variables) if variables else [None]
    dirs = None
    mode = "enter"
    if input_cmd.get("workspace", None) is not None:
        w = input_cmd.get("workspace")
        dirs = w.get("folders", None)
        if isinstance(dirs, str):
            dirs = [dirs]
        if not isinstance(dirs, list):
            raise TypeError(f"Invalid type for 'dirs': {type(dirs)}")
        mode = w.get("mode", "enter")
        qs.folder_handler(mode)
    list_of_commands = input_cmd.get("commands", None)
    if list_of_commands is None:
        raise ValueError("No commands provided in the input YAML file")
    if isinstance(list_of_commands, str):
        list_of_commands = [list_of_commands]
    elif not isinstance(list_of_commands, list):
        raise ValueError("Commands should be a list of strings or a single string")
    qs.run_script(list_of_commands, dirs, variables=list_of_variables)

if __name__ == "__main__":
    main()
