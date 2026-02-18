This documentation details the process of deploying a persistent, containerized **RStudio Server** on the Cypress HPC cluster using **Singularity** and accessing it via an **SSH Tunnel**.

---

# HPC Interactive Guide: RStudio Server via Singularity

This workflow enables a full RStudio GUI session on a compute node, utilizing a custom Conda-based R kernel and high-memory allocations.

## 1. Requesting Compute Resources

To avoid overloading the login nodes, start an interactive session. We use `idev` to request a compute node with specific port forwarding and memory.

```bash
# Requesting 250GB RAM and mapping remote port 8787 to local 18787
idev -t 10:00:00 --partition=centos7 --port=18787:8787 --mem=256000

```

---

## 2. Environment Setup

Once on the compute node, initialize your shell and activate the specialized R environment.

```bash
cd /lustre/project/sglover3/
./anaconda3/bin/conda init   # Initialize Conda for the current shell
conda activate R_env         # Activate your R environment
module load singularity/3.9.0

```

### Path Configuration

We must ensure the Singularity container can "see" the Conda binaries on the host system.

```bash
export SINGULARITYENV_PREPEND_PATH=/lustre/project/sglover3/anaconda3/bin
cd rstudio/

```

---

## 3. Launching the RStudio Container

The following command mounts necessary host directories into the container and starts the `rserver`.

```bash
singularity exec \
  -B /lustre/project/sglover3:/lustre/project/sglover3 \
  -B /lustre/project/sglover3/rstudio/tmp/var/lib:/var/lib/rstudio-server \
  -B /lustre/project/sglover3/rstudio/tmp/var/run:/var/run/rstudio-server \
  -B /lustre/project/sglover3/rstudio/tmp/tmp:/tmp \
  --workdir $(mktemp -d) \
  rstudio_latest.sif rserver \
  --www-address=0.0.0.0 \
  --server-user=$(whoami) \
  --rsession-which-r=/lustre/project/sglover3/anaconda3/envs/R_env/bin/R \
  --rsession-ld-library-path=/lustre/project/sglover3/anaconda3/envs/R_env/lib \
  --database-config-file /lustre/project/sglover3/rstudio/database.conf

```

### Key Flags Explained

| Flag | Purpose |
| --- | --- |
| `-B` (Bind) | Mounts host directories (Lustre) into the container. |
| `--rsession-which-r` | Forces RStudio to use the R version from your `R_env`. |
| `--www-address` | Listens on all interfaces within the compute node. |
| `--database-config-file` | Points to a custom SQLite database config to avoid root permission errors. |

---

## 4. SSH Tunneling & Access

To view the RStudio interface in your local browser, you must create a bridge between your computer and the HPC.

### Client-Side Configuration (PuTTY / Xshell)

1. **Host Name:** Connect to the Cypress login node as usual.
2. **Tunneling Settings:**
* **Source Port:** `8787`
* **Destination:** `localhost:18787` (or the specific compute node ID if required).


3. **Add & Open:** Click "Add" and then "Open" to start the session.

### Connecting via Web Browser

Open your preferred browser and navigate to:

> **URL:** `http://127.0.0.1:8787/`

Enter your HPC credentials when prompted. You are now running RStudio on a high-performance compute node.

---

## 5. Termination

When finished:

1. Save your workspace in RStudio.
2. Click the **Power icon (Quit)** in the top right of the RStudio UI.
3. In your terminal, use `Ctrl+C` to stop the Singularity instance.
4. Type `exit` to release the `idev` compute node allocation.

**Would you like me to help you configure the `database.conf` file mentioned in the launch command to ensure persistent user sessions?**
