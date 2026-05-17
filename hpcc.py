import os
from subprocess import call
from datetime import date

import eletility


def hpcc():
    today = date.today()
    day = today.strftime("%b-%d-%Y")

    print(
        "POET is running on the HPCC experiment mode. Please enter the following information\n"
        "======================================"
    )
    title = input("Experiment Title: ")
    hours = int(input("Requested Hours: "))
    reps = int(input("Number of Experiment Repeats: "))
    runs = int(input("Evolutionary Generations: "))
    seed = int(input("Starting Random Seed: "))
    config = input("Experiment configuration file: ")
    confirmation_text = (
        "#======================================\n#Experiment Title: {}\n#Hours: {}\n"
        "#Repeats: {}\n#Generations: {}\n#Seed: {}\n#Config: {}".format(
            title, hours, reps, runs, seed, config
        )
    )

    print(confirmation_text)

    confirm = input("To confirm the above settings, enter YES: ")
    if confirm.lower() != "yes":
        print("Aborting!")
        exit()

    title = day + "-" + title

    if os.path.exists("output/{}".format(title)):
        raise FileExistsError(
            "An experiment with the same title already exists: output/{}".format(title)
        )

    directory = os.path.join("output", title)
    subs_directory = os.path.join(directory, "subs")
    slurms_directory = os.path.join(directory, "slurms")
    logs_directory = os.path.join(directory, "logs")
    models_directory = os.path.join(directory, "models")
    errors_directory = os.path.join(directory, "errors")
    os.makedirs(directory)
    os.makedirs(subs_directory)
    os.makedirs(logs_directory)
    os.makedirs(models_directory)
    os.makedirs(errors_directory)
    os.makedirs(slurms_directory)

    content = confirmation_text
    with open(config, "r") as cfile:
        content += "\n" + cfile.read()

    file_handler = eletility.Files()
    file_handler.writeTruncate(os.path.join(directory, "config.ini"), content)

    for i in range(reps):
        filename = os.path.join(subs_directory, "{}_{}.sb".format(title, i))
        log_file_path = os.path.join(logs_directory, "evo_{}.csv".format(i))
        model_file_path = os.path.join(models_directory, "model_{}.csv".format(i))
        sbatch_content = (
            "#!/bin/bash --login\n"
            "\n########## SBATCH Lines for Resource Request ##########\n\n"
            "#SBATCH --time={}:02:00\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --mem-per-cpu=8G\n"
            "#SBATCH --job-name {}_{}\n"
            "#SBATCH --error={}/{}_{}.err\n"
            "#SBATCH --output={}/{}_{}.txt\n"
            "\n########## Command Lines to Run ##########\n\n"
            "module purge\n"
            "module load  Python/3.13.1-GCCcore-14.2.0\n"
            "cd ~/POET\n"
            "source venv/bin/activate\n"
            "srun -n 1 python poet.py -config {} -r {} -o {} -mo {} -seed {}\n"
            "cd batch\n"
            "scontrol show job $SLURM_JOB_ID\n"
        ).format(
            hours,
            title,
            i,
            errors_directory,
            title,
            i,
            slurms_directory,
            title,
            i,
            config,
            runs,
            log_file_path,
            model_file_path,
            i + seed,
        )
        with open(filename, "w") as file:
            file.write(sbatch_content)
        call(["sbatch", filename])
