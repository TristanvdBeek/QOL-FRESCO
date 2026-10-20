"""
ManagerMultiSeededMD
Tristan van der Beek
date:
This code acts as an extension onto FRESCO (by hjwijma) and pyFRESCO (by kt-korbeld) to
distribute a large number of molecular dynamics simulations for yasara for a linux PC, while keeping checkpoints
in case the program terminates. This program can help running the MD on a dedicated HPC, but was originally
written for a normal pc.
"""

DEBUG = False #can be changed by the user to get a more verbose output


import yaml #type:ignore
import pprint
import os
import sys
import shutil
import time
import datetime
import subprocess
import math
import re



def checkerror(conditional: bool, errormessage: str) -> None:
    '''
    From pyFRESCO by kt-korbeld
    Takes in conditional statement (Bool) and error message (str).
    Prints error message and quits if conditional is True
    '''
    # print message and end script if condition is met
    if conditional:
        print(f"\n-------- ERROR --------\n{errormessage}")
        sys.exit()




def prompt_user(message: str, valid_answers: tuple[str, ...]) -> str:
    """
    Prompts user with the message. Then the user can input a string
    and only valid_answers are accepted. Returns the user input
    """
    print(f"{message} ({", ".join(valid_answers)})")
    while True:
        user_input = input()
        if user_input == "":
            continue
        if user_input in valid_answers:
            return user_input
        else:
            print(f"{user_input} is not a valid answer")
            print(f"{message} ({", ".join(valid_answers)})")

        


def string_to_timedelta(time_str: str) -> datetime.timedelta:
    """
    Converts a string to a datetime.timedelta based on the conventions of this program
    """
    try:
        hours, minutes, seconds = [int(element) for element in time_str.split(":")]
        return datetime.timedelta(hours = hours, minutes = minutes, seconds = seconds) #Use format h+:mm:ss
    except:
        print(f"\n-------- ERROR --------\nThe input time {time_str} is in invalid")
        sys.exit()
        return 0




def timedelta_to_string(timedelta: datetime.timedelta) -> str:
    """
    Converts a datetime.timedelta to a string based on the conventions of this program
    """
    checkerror(type(timedelta) != datetime.timedelta, f"{timedelta} is not a datetime.timedelta")
    total_seconds = round(math.fabs(timedelta.total_seconds())) #Use format h+:mm:ss
    hours = total_seconds // 3600
    mins = (total_seconds % 3600) // 60
    seconds = round(total_seconds % 60)
    return f"{hours:02d}:{mins:02}:{seconds:02}"




def datetime_to_string(date_time):
    """
    Converts a datetime.datetime object into a string based on the conventions of this program
    """
    if date_time == None:
        return None
    elif type(date_time) == datetime.datetime:
        return date_time.strftime("%d/%m/%Y %H:%M:%S") #Use format "dd/mm/yyyy hh:mm:ss"
    else:
        checkerror(True, f"{date_time} is not a valid type for datetime_to_string conversion (type: {type(date_time)})")




def string_to_datetime(string: str | None) -> datetime.datetime | None:
    """
    Converts a string into a datetime.datetime object based on the conventions of this program
    """
    if string == None:
        return None
    elif type(string) == str:
        try:
            return datetime.datetime.strptime(string, "%d/%m/%Y %H:%M:%S") #Use format "dd/mm/yyyy hh:mm:ss"
        except:
            checkerror(True, f"{string} is not a valid format for string_to_datetime conversion")




def example_use() -> None:
    print("This is script written to automate the process of running MD simulations in YASARA "
          "after having run pyFRESCO (by kt-korbeld) or FRESCO (by hjwijma) and potentially fetchPDBs (from this repository).\n"
          "The program has the same dependencies as pyFRESCO plus the dependency on pyyaml.\n"
          "To run this program use a command in the form of: \n\tpython ManagerMultiSeededMD.py /path/to/template.pdb /path/to/mutation/pdb/directory"
          "/path/to/FRESCO/MultiSeededMD.mcr /path/to/clean_pdb.mcr <number_of_parallel_MD_runs> <max_number_failed_MD_before_give_up> <max_run_time_in hh:mm:ss> /path/to/yasara.exe \n" 
          "And example of this command is:\n\tpython ManagerMultiSeededMD.py ./KitG_dimer_AF.pdb ./QresFEP_mutations/ /software/yasara/mcr/MultiSeededMD.mcr"
          "5 3 100:04:30 /software/yasara\n"
          "The output of this program will be dumped in the directory from which this program was called.")
    return None




def parse_cli() -> dict:
    """
    Function parses the output of the ManagerMultSeededMD call. Returns a settings dictionary
    with data structured the same as what is expected of the chckpnt file.
    Function does not parse the time directly. That responsibility lies with Controller
    """
    nr_of_arguments = 8
    if DEBUG:
        print(f"The CLI command used to call this program was: \npython {" ".join(sys.argv)}")
    
    checkerror(len(sys.argv) > nr_of_arguments + 1, f"Too many arguments were passed into the function: {len(sys.argv)-1} instead of {nr_of_arguments}")
    checkerror(len(sys.argv) < nr_of_arguments + 1, f"Too few arguments were passed into the function: {len(sys.argv)-1} instead of {nr_of_arguments}")
    try:
        settings_dict = {
            "path_directory_pdb": sys.argv[2], 
            "path_template_pdb": sys.argv[1], 
            "path_multiseeded_mcr": sys.argv[3], 
            "path_clean_pdb_mcr": sys.argv[4],
            "path_yasara": sys.argv[8], 
            "number_of_parallel_runs": int(sys.argv[5]), 
            "max_number_of_reruns_failed_md": int(sys.argv[6]), 
            "max_run_time": sys.argv[7]
        }
    except:
        print(f"\n-------- ERROR --------\npython {" ".join(sys.argv)} is not a valid command")
        sys.exit()

    checkerror(not os.path.isdir(settings_dict["path_directory_pdb"]), f'"{settings_dict["path_directory_pdb"]} is not a valid directory')
    checkerror(not os.path.isfile(settings_dict["path_template_pdb"]) and settings_dict["path_template_pdb"][-4] == ".pdb", f'"{settings_dict["path_template_pdb"]} is not a valid pdb file')
    checkerror(not os.path.isfile(settings_dict["path_multiseeded_mcr"]) and settings_dict["path_multiseeded_mcr"][-4] == ".mcr", f'"{settings_dict["path_multiseeded_mcr"]} is not a valid macro file')
    checkerror(not os.path.isfile(settings_dict["path_clean_pdb_mcr"]) and settings_dict["path_clean_pdb_mcr"][-4] == ".mcr", f'"{settings_dict["path_clean_pdb_mcr"]} is not a valid macro file')
    #FIGURE OUT HOW TO DO THIS LINE
    checkerror(os.access(settings_dict["path_yasara"], os.X_OK) == False, f"{settings_dict["path_yasara"]} is not a valid path to an executable file")
    # checkerror(os.path.isfile(settings_dict["path_yasara"]), f'"{settings_dict["path_template_pdb"]} is not a valid pdb file')
    checkerror(len(settings_dict["max_run_time"].split(":")) != 3 or any([True for time in settings_dict["max_run_time"].split(":")[1:] if len(time) != 2]),
                f"{settings_dict["max_run_time"]} is not a valid format for time")
    if DEBUG:
        print("The following settings were passed to the program")
        pprint.pp(settings_dict)
    return settings_dict




class Chckpnt_parser:
    _path_checkpoint_file = os.path.join(".", "chckpnt.yaml")
    _path_temp_checkpoint_file = _path_checkpoint_file[:-5] + "_temp_" + _path_checkpoint_file[-5:]
    _entries_yaml = ("settings", "metadata", "previous_settings", "data_md_runs")
    _entries_settings = ("path_directory_pdb", "path_template_pdb", "path_multiseeded_mcr", "path_clean_pdb_mcr", "path_yasara", 
                    "number_of_parallel_runs", "max_number_of_reruns_failed_md", "max_run_time")
    _entries_metadata = ("start_datetime", "end_datetime", "run_time", "status")

    def _safe_write(self, chckpnt: dict) -> None:
        """
        Used to dump a chckpnt to the yaml file in a secure manner, that should not corrupt the file when an interrupt
        or shutdown occurs
        """
        with open(self._path_temp_checkpoint_file, "w") as open_file:
            yaml.dump(chckpnt, open_file, sort_keys=False)
            open_file.flush()
            os.fsync(open_file.fileno())
        os.replace(self._path_temp_checkpoint_file, self._path_checkpoint_file)


    def read_previous_chckpnt(self) -> None | dict:
        """
        Reads the checkpoint file, returns None if it doesn't exist,
        if it does exist, the function checks whether all required fields are there.
        If all fields are there as expected, the function outputs a dictionary with data 
        from the \"settings\" entry, if fields are missing, the program is
        terminated.
        """
        if os.path.isfile(self._path_checkpoint_file) == False:
            if DEBUG:
                print(f"No checkpoint file was found at {self._path_checkpoint_file}")
            return None
        
        #This gets executed if the checkpoint file exists
        with open(self._path_checkpoint_file, "r") as open_file:
            chckpnt = yaml.safe_load(open_file)
        if DEBUG:
            print(f"Found a checkpoint file at {self._path_checkpoint_file}")

        chckpnt_checks = {
            "general categories": not False in [entry in list(chckpnt.keys()) for entry in self._entries_yaml],
            "settings":  not False in [entry in list(chckpnt["settings"].keys()) for entry in self._entries_settings],
            "metadata": not False in [entry in list(chckpnt["metadata"].keys()) for entry in self._entries_metadata]
        }
        
        if DEBUG:
            print("The found checkpoint file performed as followed on the following checks", 
                  f"\n{"\n".join([str(key) + ": " + str(value) for key, value in chckpnt_checks.items()])}")

        checkerror(False in tuple(chckpnt_checks.values()), f"{self._path_checkpoint_file} in the directory is invalid")
        return chckpnt
    

    def read_chckpnt(self) -> dict:
        """
        Reads and outputs the entire chckpnt
        """
        with open(self._path_checkpoint_file, "r") as open_file:
            chckpnt = yaml.safe_load(open_file)
        if DEBUG:
            print(f"The manager data from the chckpnt file was read")
        return chckpnt


    def read_md(self) -> dict:
        """
        Read the chckpnt file, outputs a dictionary with the data from the 
        data_md_runs section
        """
        if DEBUG:
            print(f"The MD data was read from the chckpnt file")
        return self.read_chckpnt()["data_md_runs"]
    

    def write_to_chckpnt(self, chckpnt: dict) -> None:
        """
        Writes a chckpnt to the chckpnt file, also updates the end_time
        """
        if DEBUG:
            print("A complete chckpoint will be written")

        chckpnt["metadata"]["end_datetime"] = datetime_to_string(datetime.datetime.now())
        chckpnt["metadata"]["run_time"] = timedelta_to_string(string_to_datetime(chckpnt["metadata"]["end_datetime"]) - 
                                                              string_to_datetime(chckpnt["metadata"]["start_datetime"])) # type: ignore
        self._safe_write(chckpnt)


    def write_md_run(self, md_run_dict: dict) -> None: 
        """
        Writes the data in md_run_dict to the chckpnt file. Also updates end_time and run_time
        of (manager) metadata section
        The md_run_dict contains the keys: mutation, status, start_time, end_time, run_time
        """
        #Update the metadata of the mutation
        with open(self._path_checkpoint_file, "r") as open_file:
            chckpnt = yaml.safe_load(open_file)
        chckpnt["data_md_runs"][md_run_dict["mutation"]] = {
            "start_time": md_run_dict["start_time"],
            "end_time": md_run_dict["end_time"],
            "run_time": md_run_dict["run_time"],
            "status": md_run_dict["status"],
            "number_of_runs": md_run_dict["number_of_runs"]
        }

        chckpnt["metadata"]["end_datetime"] = datetime_to_string(datetime.datetime.now())
        chckpnt["metadata"]["run_time"] = timedelta_to_string(string_to_datetime(chckpnt["metadata"]["end_datetime"]) - 
                                                              string_to_datetime(chckpnt["metadata"]["start_datetime"])) #type:ignore
        
        if DEBUG:
            print(f"{md_run_dict["mutation"]} has its data dumped into the chckpnt file, chckpnt metadata was updated:")
            pprint.pp(chckpnt["data_md_runs"][md_run_dict["mutation"]])
        self._safe_write(chckpnt)
        return None

    def write_updated_end_time(self) -> None:
        with open(self._path_checkpoint_file, "r") as open_file:
                    chckpnt = yaml.safe_load(open_file)
        
        chckpnt["metadata"]["end_datetime"] = datetime_to_string(datetime.datetime.now())
        chckpnt["metadata"]["run_time"] = timedelta_to_string(string_to_datetime(chckpnt["metadata"]["end_datetime"]) - 
                                                                string_to_datetime(chckpnt["metadata"]["start_datetime"])) #type:ignore
        self._safe_write(chckpnt)
        return None
        




class Yasara_runner:
    PATH_OUTPUT_DIR = "./"
    TIME_OUT = datetime.timedelta(hours = 1, minutes = 30)
    def __init__(self, path_yasara, path_multiseeded_md, path_clean_pdb: str):
        #these variables are not changed after initiation
        self.path_yasara = path_yasara
        self.path_multiseeded_md = path_multiseeded_md
        self.path_clean_pdb = path_clean_pdb

    def new_md(self, path_structure: str, mutation: str, number_of_runs: int) -> None:
        if number_of_runs > 0:
            print(f"Going to run MD number {number_of_runs} for {mutation}")
        else:
            print(f"Going to run the MD for {mutation}")

        self._initialise_run(path_structure, mutation, number_of_runs)
        self._prepare_directory()
        self.command = self._prepare_command()
        self._start_process(self.command)


    def _initialise_run(self, path_structure: str, mutation: str, number_of_runs: int) -> None:
        self.start_time = datetime.datetime.now()
        self.end_time = None
        self.disabled = False
        self.path_structure_old_dir = path_structure
        self.mutation = mutation
        self.number_of_runs = number_of_runs
        self.path_dir = os.path.join(self.PATH_OUTPUT_DIR, f"Subdir_{mutation}")
        self.name_structure_for_mcr = self._write_name_mcr()
        self.path_structure = os.path.join(self.path_dir, self.name_structure_for_mcr + ".pdb")


    def _write_name_mcr(self):
        name = os.path.split(self.path_structure_old_dir)[-1]
        if self.mutation == "template":
            prefix = name.split(".")[0]
            print("prefix", prefix) #REMOVE
            return f"{prefix}_cleaned_{self.mutation}"
        else:
            prefix = name.split(self.mutation)[0]
            return f"{prefix}cleaned_{self.mutation}_WOW"


    def _prepare_directory(self) -> None:
        """
        Empties self.path_dir or makes the directory
        Then add the pdb file of the mutation to the directory        
        """
        dir_exists = os.path.isdir(self.path_dir)

        #Remove directory and its contents 
        if dir_exists:
            files_in_dir = os.listdir(self.path_dir)
            if len(files_in_dir) > 0: 
                print(f"{self.path_dir} already exists, so its contents will be removed")
                if DEBUG:
                    print(f"{self.path_dir} contained: {", ".join(files_in_dir)}")
            for file_name in os.listdir(self.path_dir):
                if file_name[:3] == "LOG" or file_name[:-4] == self.name_structure_for_mcr:
                    pass
                else:
                    os.remove(os.path.join(self.path_dir, file_name))
        else:
        #Making the directory and copying the pdb file into it
            os.mkdir(self.path_dir)
            shutil.copyfile(self.path_structure_old_dir, self.path_structure)
            self._clean_structure()


    def _clean_structure(self) -> None:
        if DEBUG:
            print(f"Cleaning up {self.name_structure_for_mcr}")
        command = f"cd {self.path_dir}; {self.path_yasara} -txt {self.path_clean_pdb} \"MacroTarget = \'{self.name_structure_for_mcr}\'\" > LOG_CLEAN_PDB"
        subprocess.run(command, shell=True)
        return


    def _prepare_command(self) -> str:
        command =  f"cd {self.path_dir};{self.path_yasara} -txt {self.path_multiseeded_md} \"MacroTarget = \'{self.name_structure_for_mcr}\'\" > LOG"
        if DEBUG: 
            print(f"The following command is now being run: {command}")
        return command

    def _start_process(self, command) -> None:
        self.process = subprocess.Popen(command, shell = True)


    def check_has_finished(self) -> int:
        #Needs to return 0 or 1 depending on if it has finished or not, -1 is failed, 2 is yasara_runner is disabled
        return_code = self.process.poll()
        self._check_overtime()
        if DEBUG:
            print(f"{self.mutation} gave return code: {return_code}")
        if self.disabled == True:
            return 2
        if return_code == None:
            return 0
        elif return_code >= 0 and self._check_succesful_md_run() == True:
            if self.end_time == None:
                self.end_time = datetime.datetime.now()
                self.number_of_runs += 1
                if self.number_of_runs > 1: 
                    print(f"{self.mutation} finished running succesfully after having failed {self.number_of_runs-1} runs")
                else:
                    print(f"{self.mutation} finished running succesfully")
                if DEBUG:
                    print(f"The MD simulation took {timedelta_to_string(self.end_time - self.start_time)}")
            return 1
        elif return_code < 0 or self._check_succesful_md_run() == False:
            if self.end_time == None:
                self.end_time = datetime.datetime.now()
                self.number_of_runs += 1
                print(f"{self.mutation} failed its MD simulation number {self.number_of_runs}")
                if DEBUG:
                    print(f"The MD simulation took {timedelta_to_string(self.end_time - self.start_time)}")
            return -1
        else:
            print(f"It just outputs return code {self.mutation}")
            return return_code 


    def _check_overtime(self) -> None:
        if (datetime.datetime.now() - self.start_time > Yasara_runner.TIME_OUT) == True:
            self.process.kill()
            print(f"{self.mutation} ran overtime and was thus killed")
        return 


    def get_metadata(self) -> dict:
        """
        Function used to get metadata of a running or finished MD simulation
        """
        #Used to get data such as mutation name, time elapsed, etc. failed (interrupt should be handled by controller)
        return_code = self.check_has_finished()
        status: str = ""
        if return_code == 0:
            status = "running"
        elif return_code == 1:
            status = "finished"
        elif return_code == -1:
            status = "failed (program_failure)"
        else:
            checkerror(True, f"Wrong return code: {return_code}")

        metadata_dict = {
            "mutation": self.mutation,
            "status": status,
            "start_time": datetime_to_string(self.start_time),
            "end_time": datetime_to_string(self.end_time),
            "run_time": None,
            "number_of_runs": self.number_of_runs
        }

        if self.end_time != None:
            metadata_dict["run_time"] = timedelta_to_string(self.end_time - self.start_time)
        return metadata_dict                


    def _check_succesful_md_run(self) -> bool:
        """
        Checks whether the MD was successful. Return True if <mutation>_All5Seeds_LSOn_AvgRMSF.pdb was there
        and False if it is absent.
        """
        #Returns True succesful, False if unsuccesful
        # n 1 #CHANGE, UNCOMMENT
        ENDING_IMPORTANT_FILE = "All5Seeds_LSOn_AvgRMSF.pdb" #CHECK
        files_dir = os.listdir(self.path_dir)
        for file in files_dir:
            if ENDING_IMPORTANT_FILE in file:
                return True
        return False


    def rerun_md(self) -> None:
        """
        Reruns the MD simulation with the exact same parameters as were given for the previous run
        """
        print(f"Going to rerun the MD simulation for mutation {self.mutation}") #CHANGE
        shutil.copyfile(os.path.join(self.path_dir, "LOG"), os.path.join(self.path_dir, f"LOG_run_{self.number_of_runs}"))

        self._prepare_directory()
        self._initialise_run(self.path_structure_old_dir, self.mutation, self.number_of_runs)
        self._start_process(self.command)


    def disable(self) -> None:
        """
        Disables the Yasara_runner. If this function is called, all check_has_finished will return 2 and
        get_metadata can no longer be called.
        """
        self.disabled = True


      

class Controller:
    TIME_BETWEEN_LOADING_BARS = 150 #seconds
    TIME_BETWEEN_CHECKS = 10 #seconds
    def __init__(self, settings_dict: dict) -> None:
        """
        Master controller of the manager. This calls all the classes, functions and methods (except for Cli_parser) to manage the simulation
        """
        #Setup a new chckpnt file with the correct settings, data and initialised md_runs. If a chckpnt already exists, updates the chckpnt file instead
        self.settings_dict: dict = settings_dict
        self.md_run_dict = {}
        self.previous_settings: list = []
        self.chckpnt_parser = Chckpnt_parser()
        self.previous_chckpnt: dict | None = self.chckpnt_parser.read_previous_chckpnt()
        if self.previous_chckpnt is not None:
            self.check_match_settings()
        else:
            self.make_new_settings()

        self.start_time = datetime.datetime.now()
        self.number_of_disabled_yasara_runners = 0
        
        self.metadata = {
            "start_datetime": datetime_to_string(self.start_time),
            "end_datetime": datetime_to_string(datetime.datetime.now()),
            "run_time": timedelta_to_string(datetime.datetime.now() - self.start_time),
            "status": "running"
        }
        
        chckpnt = {
            "settings": self.settings_dict,
            "metadata": self.metadata,
            "previous_settings": self.past_runs_metadata,
            "data_md_runs": self.md_run_dict
        }

        self.chckpnt_parser.write_to_chckpnt(chckpnt)

        #Initialise the yasara runners and active all of them
        self.yasara_runners: list[Yasara_runner] = [Yasara_runner(self.settings_dict["path_yasara"], self.settings_dict["path_multiseeded_mcr"], self.settings_dict["path_clean_pdb_mcr"]) for i in range(settings_dict["number_of_parallel_runs"])]
        #For every yasara runner, start a new fresh run and update the chckpnt file
        for yasara_runner in self.yasara_runners:
            self.find_new_md(yasara_runner)
            self.update_chckpnt_md(yasara_runner)

        time_previous_loading_bar = datetime.datetime.now()

        #Enter the main managing loop
        while True: #update
            #Check if passed max_run_time and check if all yasara runners are done
            self.chckpnt_parser.write_updated_end_time()
            if self.number_of_disabled_yasara_runners == self.settings_dict["number_of_parallel_runs"]:
                self.finish_up_succesful()
                break
            elif self.time_left().total_seconds() < 0:
                self.finish_up_overtime()
                break
            #Loop over yasara runners and check their status. If needed activate a new MD or rerun an old md
            for yasara_runner in self.yasara_runners:
                if yasara_runner.check_has_finished() in [0,2]:
                    continue
                #If the run failed and it hasn't reached the maximum number of runs yet, perform a rerun
                elif yasara_runner.check_has_finished() == -1 and yasara_runner.get_metadata()["number_of_runs"] < self.settings_dict["max_number_of_reruns_failed_md"]:
                    yasara_runner.rerun_md()
                    self.update_chckpnt_md(yasara_runner)
                    continue
                #do a new mutation
                elif yasara_runner.check_has_finished() == 1 or (yasara_runner.get_metadata()["number_of_runs"] >= self.settings_dict["max_number_of_reruns_failed_md"]):
                    self.update_chckpnt_md(yasara_runner)
                    self.find_new_md(yasara_runner)
                    if yasara_runner.check_has_finished() != 2: #To prevent asking for metadata when the yasara_runner has been disabled
                        self.update_chckpnt_md(yasara_runner)
            #If enough time has elapsed, write a new loading bar
            if (datetime.datetime.now() - time_previous_loading_bar).total_seconds() > self.TIME_BETWEEN_LOADING_BARS:
                self.write_loading_bar()
                time_previous_loading_bar = datetime.datetime.now()

            time.sleep(self.TIME_BETWEEN_CHECKS)
        

    def update_chckpnt_md(self, yasara_runner: Yasara_runner):
        """
        Update the chckpnt with new data of the input Yasara_runner
        """
        metadata = yasara_runner.get_metadata()
        self.chckpnt_parser.write_md_run(metadata)
        return None

            
    def find_new_md(self, yasara_runner: Yasara_runner) -> None:
        """
        Finds a new mutation for MD in the chckpnt file, then intialises this run in the
        input Yasara_runner and update the chckpnt file.
        If there are no more mutations to run, disable the yasara runner instead
        """
        md_runs = self.chckpnt_parser.read_md()
        #Find the next mutation that has not yet been run
        for mutation, metadata in md_runs.items(): 
            if metadata == None:
                self._update_yasara_runner(yasara_runner, mutation, 0)
                break
            elif "failed" in metadata["status"] and metadata["number_of_runs"] < self.settings_dict["max_number_of_reruns_failed_md"]:
                self._update_yasara_runner(yasara_runner, mutation, metadata["number_of_runs"])
                break  
        else: #If no new mutations are found, disable the yasara runner 
            yasara_runner.disable()
            self.number_of_disabled_yasara_runners += 1
            return None
        return None


    def _update_yasara_runner(self, yasara_runner: Yasara_runner, mutation: str, number_of_runs: int) -> None:
        if mutation == "template":
            path_pdb: str = self.settings_dict["path_template_pdb"]
        else:
            path_pdb: str = [os.path.join(self.settings_dict["path_directory_pdb"], file) for file in os.listdir(self.settings_dict["path_directory_pdb"]) if file[-len(mutation)-4:-4] == mutation][0]
        yasara_runner.new_md(path_pdb, mutation, number_of_runs)
        return None



    def check_match_settings(self) -> None:
        """
        This function is run if there already exists a chckpnt file. 
        It checks whether the parameters max_number_of_reruns_failed_md, path_directory_pdb, path_template_pdb, path_multiseeded_mcr, path_clean_pdb_mcr are different
        if they, prompt the user. Then it checks if the mutations in the pdb directory and in the chckpnt file match. lastly, updates the status of
        the mutations from running to failed (interrupt), the status of the previous program to failed (interrupt) if it was at running. Then puts 
        the metadata of the program into past_runs_metadata
        """
        assert(self.previous_chckpnt is not None)
        #Checking whether and prompting if the max number of reruns failed md has changed
        if self.settings_dict["max_number_of_reruns_failed_md"] != int(self.previous_chckpnt["settings"]["max_number_of_reruns_failed_md"]):
            if prompt_user("You changed the maximum number of reruns after a failed MD simulation before the program stops rerunning MD "
                        + f"from {self.previous_chckpnt["settings"]["max_number_of_reruns_failed_md"]} to {str(self.settings_dict["max_number_of_reruns_failed_md"])}\n"
                        + "do you want to override all previous MD simulations with the new parameters?", ("y", "n")) == "n":
                sys.exit()
        #Checking whether and prompting if the path_directory_pdb, path_template_pdb, path_multiseeded_mcr, path_clean_pdb_mcr failed have changed
        match_settings_dict_non_critical = {
            "path_directory_pdb": self.settings_dict["path_directory_pdb"] != self.previous_chckpnt["settings"]["path_directory_pdb"],
            "path_template_pdb": self.settings_dict["path_template_pdb"] != self.previous_chckpnt["settings"]["path_template_pdb"],
            "path_multiseeded_mcr": self.settings_dict["path_multiseeded_mcr"] != self.previous_chckpnt["settings"]["path_multiseeded_mcr"],
            "path_clean_pdb_mcr": self.settings_dict["path_clean_pdb_mcr"] != self.previous_chckpnt["settings"]["path_clean_pdb_mcr"]
        }
        
        if any(match_settings_dict_non_critical.values()):
            if prompt_user(f"The following settings were changed {", ".join([key for key, value in match_settings_dict_non_critical.items() if value == True])} " 
                        + "do you still want to continue running MD simulations?", ('y', 'n')) == "n":
                sys.exit()
        
        #Go over the mutations in the directory and in the chckpnt file and see if they match, other wise prompt user and update the chckpnt ifle
        pdb_files = [file for file in os.listdir(self.settings_dict["path_directory_pdb"]) if file[-4:].lower() == ".pdb"]
        mutations_chckpnt = tuple(self.previous_chckpnt["data_md_runs"].keys())
        if DEBUG:
            print(f"The following mutations were found in the chckpnt file:\n{", ".join(mutations_chckpnt)}")
        mutations_for_from_files = [mutation.split("_")[-1][:-4] for mutation in pdb_files] + ["template"]
        new_mutations = [mutation for mutation in mutations_for_from_files if mutation not in mutations_chckpnt]
        if DEBUG:
            print(f"The following mutations were found in {self.settings_dict["path_directory_pdb"]}:\n{", ".join(mutations_for_from_files)}")
        old_mutations_not_in_directory = [mutation for mutation in mutations_chckpnt if mutation not in mutations_for_from_files]

        #Prompt user and update chckpnt file
        if len(new_mutations) > 0 or len(old_mutations_not_in_directory) > 0:
            if prompt_user(f"There is a mismatch between the mutations in the chckpnt file and the mutations in the pdb directory.\n"
                        + f"{", ".join(new_mutations)} are not in the chckpnt file and {", ".join(old_mutations_not_in_directory)} are not in the directory.\nContinuing means "
                        + f"that you gain the mutations in the directory and lose the mutations that are not in the directory. Do you want to continue", ("y", "n")) == "n":
                sys.exit()
            else: 
                self.md_run_dict: dict = self.previous_chckpnt["data_md_runs"]
                #Add new mutations from the pdb directory
                for mutation in new_mutations:
                    self.md_run_dict[mutation] = None
                #Remove mutations not in the pdb directory
                for mutation in old_mutations_not_in_directory:
                    self.md_run_dict.pop(mutation)
        else:
            self.md_run_dict: dict = self.previous_chckpnt["data_md_runs"]

        #Changes the status of all the mutations with the "running" status to failed <number_of_error_attempts> (interrupt)
        for mutation, metadata in self.md_run_dict.items(): #type:ignore
            if type(metadata) == dict and "status" in metadata.keys() and metadata["status"] == "running":
                self.md_run_dict[mutation]["status"] = "failed (interrupt)"
                self.md_run_dict[mutation]["number_of_runs"]

        #Update status of the Manager if it was running to failed (interrupt)
        status_previous = self.previous_chckpnt["metadata"]["status"]
        checkerror(self.previous_chckpnt["metadata"]["status"] == "finished" and (len(new_mutations) > 0 or len(old_mutations_not_in_directory) > 0), "The MD simulations have already been finished")
        #If the previous settings still show "running" it will replace it with failed (interrupt) in the previous_metadata. Thus any finished calculation needs to set finished or failed (timeout) itself
        if self.previous_chckpnt["metadata"]["status"] == "running":
            self.previous_chckpnt["metadata"]["status"] = "failed (interrupt)"
        
        #Push previous settings onto self.past_runs_metadata 
        self.past_runs_metadata = self.previous_chckpnt["previous_settings"]
        new_previous_metadata = {
            "settings": self.previous_chckpnt["settings"],
            "metadata": self.previous_chckpnt["metadata"]
        }
        self.past_runs_metadata.insert(0,new_previous_metadata)
        return None


    def make_new_settings(self) -> None:
        """
        If no chckpnt file, initialise some important data
        """
        mutations = [mutation.split("_")[-1][:-4] for mutation in os.listdir(self.settings_dict["path_directory_pdb"]) if mutation[-4:].lower() == ".pdb"] + ["template"]
        print(mutations) #remove
        self.md_run_dict = {mutation: None for mutation in mutations}
        self.past_runs_metadata = []
        return None


    def write_loading_bar(self):
        """
        Print a loading bar to the terminal based on how many MD runs have been completed and how much time is remaining
        """
        chckpnt = self.chckpnt_parser.read_chckpnt()
        nr_mds = len(chckpnt["data_md_runs"].keys())
        nr_finished_mds = len([mutation for mutation, metadata in chckpnt["data_md_runs"].items() if metadata != None and (metadata["status"] == "finished" or metadata["number_of_runs"] >= self.settings_dict["max_number_of_reruns_failed_md"])])
        percentage_finished = round(nr_finished_mds / nr_mds * 100)
        time_remaining: str = timedelta_to_string(self.time_left())
        print(f"\n{time_remaining} out of {chckpnt["settings"]["max_run_time"]} time remaining\n" +
              f"[{"#" * (percentage_finished // 4)}{"." * (25 - percentage_finished // 4)}] {str(percentage_finished)} %", flush = True)
            

    def time_left(self) -> datetime.timedelta:
        """
        Returns a datetime.timedelta with the amount of time left for the program
        """
        time_elapsed: datetime.timedelta = datetime.datetime.now() - string_to_datetime(self.metadata["start_datetime"]) #type:ignore
        return string_to_timedelta(self.settings_dict["max_run_time"]) - time_elapsed


    def finish_up_overtime(self):
        """
        Used to finish up writing data to the chckpnt file in case of overtime
        """
        chckpnt = self.chckpnt_parser.read_chckpnt()
        chckpnt["metadata"]["end_datetime"] = datetime_to_string(datetime.datetime.now())
        chckpnt["metadata"]["run_time"] = timedelta_to_string(string_to_datetime(chckpnt["metadata"]["end_datetime"]) - 
                                                              string_to_datetime(chckpnt["metadata"]["start_datetime"])) #type:ignore
        chckpnt["metadata"]["status"] = "failed (overtime)"
        self.chckpnt_parser.write_to_chckpnt(chckpnt)

        md_data = self.chckpnt_parser.read_md()
        nr_mds = len(md_data.keys())
        failed_mds = [mutation for mutation, metadata in md_data.items() if metadata != None and metadata["status"] == "failed (program_failure)"]
        unfinished_mds = [mutation for mutation, metadata in md_data.items() if metadata == None or metadata["status"] not in ["failed (program_failure)", "finished"]]
        
        nr_failed_mds = len(failed_mds)
        print("\n-------- OVERTIME ERROR --------\n" +
            "The program finished due to running longer than the indicated max_run_time.\n" + 
            f"The program ran for {chckpnt["metadata"]["run_time"]}, while the max_run_time was {chckpnt["settings"]["max_run_time"]}. \n"
            f"Molecular dynamics was run for {nr_mds - len(unfinished_mds)} mutations of which {nr_failed_mds} failed more than {self.settings_dict["max_number_of_reruns_failed_md"]} times. "
            f"{len(unfinished_mds)} mutations are not yet finished running. The failed mutations were:\n" + ", ".join(failed_mds))


    def finish_up_succesful(self):
        """
        Used to finish up writing data to the chckpnt file in case of all simulations being finished
        """
        chckpnt = self.chckpnt_parser.read_chckpnt()
        time_run = chckpnt["metadata"]["run_time"]
        chckpnt["metadata"]["status"] = "finished"
        self.chckpnt_parser.write_to_chckpnt(chckpnt)

        md_data = self.chckpnt_parser.read_md()
        nr_mds = len(md_data.keys())
        failed_mds = [mutation for mutation, metadata in md_data.items() if metadata["status"] == "failed (program_failure)"]
        nr_failed_mds = len(failed_mds)
        print("\n-------- END OF PROGRAM --------\n" +
            "Finished running all molecular dynamics simulations\n" + 
            f"The program ran for {time_run}"
            f"Molecular dynamics was run for {nr_mds} mutations of which {nr_failed_mds} failled more than {self.settings_dict["max_number_of_reruns_failed_md"]} times. These were:\n" + 
            ", ".join(failed_mds))






def main():
    example_use()
    Controller(parse_cli())


if __name__ == "__main__":
    main()
