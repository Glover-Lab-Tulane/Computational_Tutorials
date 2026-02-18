# HPC Environment Setup: Anaconda & R on Cypress

This guide outlines the process for installing Anaconda3 in a non-standard directory (Lustre parallel filesystem) and initializing a reproducible R environment via a YAML configuration.

## 1. Installation & Directory Configuration

On HPC systems like Cypress, installing large software suites in your home directory often leads to quota issues. We bypass this by installing directly to the project storage.

### Download the Installer

Retrieve the specific Linux x86_64 distribution of Anaconda3 (2025.06) using `wget`:

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2025.06-0-Linux-x86_64.sh

```

### Custom Installation Path

Execute the shell script. When prompted for the installation location, override the default path with your project directory:

* **Target Directory:** `/lustre/project/sglover3/anaconda3`

```bash
bash Anaconda3-2025.06-0-Linux-x86_64.sh

```

---

## 2. Shell Initialization

Since the installation is in a non-standard path, you must manually point to the binary to initialize your shell (e.g., bash or zsh). This command appends the necessary conda initialization code to your `~/.bashrc`.

```bash
./anaconda3/bin/conda init

```

> [!IMPORTANT]
> After running `init`, you must **restart your shell session** or run `source ~/.bashrc` for the changes to take effect and for the `conda` command to be recognized globally.

---

## 3. Environment Provisioning (R_env)

To ensure reproducibility, we use a pre-defined YAML configuration file ([`R_env.yml`](R_env.yml)). This file contains the environment name, channels, and specific R dependencies.

### Create Environment from YAML

Run the following command from the directory containing your `.yml` file:

```bash
conda env create -f R_env.yml

```

### Environment Workflow

Once the solver finishes installing the packages, manage the environment using the standard CLI:

| Action | Command |
| --- | --- |
| **Activate** | `conda activate R_env` |
| **Deactivate** | `conda deactivate` |
| **Verify Packages** | `conda list` |
