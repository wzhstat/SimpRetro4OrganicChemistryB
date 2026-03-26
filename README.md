# SimpRetro4OrganicChemistryB
A single-step retrosynthesis algorithm specifically tailored for teaching purposes.

## Installation

**1. Clone the repository**
```bash
git clone https://github.com/yourusername/SimpRetro4OrganicChemistryB.git
cd SimpRetro4OrganicChemistryB
```

**2. Create a virtual environment (Recommended)**
It is highly recommended to use `conda` to manage dependencies, especially for RDKit.
```bash
conda create -n retro_env python=3.8
conda activate retro_env
```

**3. Install dependencies**
Install RDKit, `rdchiral`, and other required packages:
```bash
pip install -r requirements.txt
```

**4. Data Preparation**
Ensure the following necessary data files are placed in the root directory of the project (or specify their paths during inference):
- `reaction_template.json`
- `template_condition.json`
- `<database>.txt` (e.g., `emol_under_0_carbons.txt`)

---

## Inference

Run the python script from your terminal. The basic usage only requires a target SMILES string.

### Basic Run
```bash
python main.py -s "CC(=O)C=C(C)C"
```
*The result will be automatically saved to `retro_result.json`.*

### Advanced Run (Customizing Weights & Templates)
You can customize the scoring weights ($w_1, w_2, w_3, w_4$) and specify your own template/condition files:
```bash
python main.py \
    -s "CC(=O)C=C(C)C" \
    -w 0.2 0.3 0.5 0.0 \
    -tpl "custom_templates.json" \
    -cond "custom_conditions.json" \
    -o "my_output.json"
```

### Command Line Arguments
| Argument | Short | Default | Description |
| :--- | :---: | :--- | :--- |
| `--smiles` | `-s` | `CC(=O)C=C(C)C` | The target molecule's SMILES string. |
| `--weights` | `-w` | `0.1 0.2 0.5 0.0` | 4 float numbers for the scoring function weights. |
| `--template` | `-tpl` | `reaction_template.json` | Path to the reaction templates JSON file. |
| `--condition`| `-cond`| `template_condition.json`| Path to the reaction conditions JSON file. |
| `--database` | `-db` | `emol_under_0` | In-stock molecules database prefix name. |
| `--output` | `-o` | `retro_result.json` | Path to save the final JSON output. |

---
*Note: Please replace `main.py` with your actual Python script filename.*
