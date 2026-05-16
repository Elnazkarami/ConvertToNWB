# ConvertToNWB

Convert Peyrache Lab postsubiculum recordings to [NWB](https://www.nwb.org/) (Neurodata Without Borders) format for publication to the DANDI archive ([dandiset 000939](https://dandiarchive.org/dandiset/000939)).

## Associated publication

Duszkiewicz AJ, *et al.* *Nature Neuroscience* (2024). [doi:10.1038/s41593-024-01588-5](https://doi.org/10.1038/s41593-024-01588-5)

## What the script does

[`PeyracheToNWB.py`](PeyracheToNWB.py) walks each session listed in `Dataset_metadata.xlsx`, builds one `.nwb` file per session, runs `nwbinspector` against it, then hands the output to `dandi organize` so it can be uploaded to the dandiset.

Each NWB file includes:

| Section | Contents |
| --- | --- |
| `behavior` | (x, y) position and head direction from Optitrack Motive 2.0; 3-axis accelerometer from the Intan RHD headstage (when present) |
| Epochs | Experimenter-defined epochs (open field, sleep, etc.) plus REM / NREM intervals from `SleepState.states.mat` |
| `ecephys` | LFP at 1250 Hz, optional raw broadband at 20 kHz, pseudo-EMG derived from LFP correlations, electrode table with faulty-channel flags |
| `units` | Spike times, mean waveforms, trough-to-peak, cell-type tags (head-direction tuned / putative excitatory / putative fast-spiking) |
| `ogen` | Optogenetic stimulus site and digital-input trace (when present) |
| Subject | Subject ID, species, mouse line, virus, age, sex |

## Expected directory layout

The script reads from `DATAPATH` (set at the top of the script). For each session the following structure is expected:

```
<DATAPATH>/
├── Dataset_metadata.xlsx        # one row per session
├── NWB/000939/                  # dandiset target directory
└── <session-folder>/            # e.g. A1234-220301-...
    ├── ChannelMap.mat           # chanOrder, goodChans
    ├── <session-folder>.dat     # raw broadband (optional, see ADD_RAW)
    ├── <session-folder>.lfp     # LFP
    ├── <session-folder>_auxiliary.dat   # accelerometer (optional)
    ├── Opto_digitalin.dat       # optogenetic digital input (optional)
    ├── Data/
    │   ├── SpikeData.mat
    │   ├── Angle.mat
    │   ├── Position.mat
    │   ├── Epoch_TS.csv
    │   ├── Waveforms.mat
    │   ├── WaveformFeatures.mat
    │   └── CellTypes.mat
    └── Sleep/
        ├── <session-folder>.SleepState.states.mat
        └── <session-folder>.EMGFromLFP.LFP.mat
```

Session folder names are expected to follow `<subject_id>-YYMMDD-...`; the script parses the date out of the second `-`-separated token.

## Metadata spreadsheet

`Dataset_metadata.xlsx` is the source of truth for per-session metadata. Required columns:

| Column | Purpose |
| --- | --- |
| `Recording` | Session folder name |
| `Virus` | Virus description (free text) |
| `Mouse_line` | Mouse line / strain |
| `Probe_description` | Free-text description of the probe |
| `Epoch_1` … `Epoch_4` | Tag string for each epoch (leave blank if absent) |
| `Opto_device`, `Opto_site_name`, `Opto_site_description`, `Opto_site_location`, `Wavelength` | Optogenetics fields (used only if `Opto_digitalin.dat` is present) |

Optional columns (the script falls back to sensible defaults if any are missing):

| Column | Default |
| --- | --- |
| `Age` | `P12W` |
| `Sex` | `M` |
| `N_shanks` | `1` |
| `Raw_sample_rate` | `20000` |
| `LFP_sample_rate` | `1250` |

## Requirements

- Python 3.9+
- [`pynwb`](https://pynwb.readthedocs.io/) and [`hdmf`](https://hdmf.readthedocs.io/)
- [`pynapple`](https://pynapple-org.github.io/pynapple/)
- `pandas`, `numpy`, `scipy`, `pytz`, `openpyxl` (for reading `.xlsx`)
- [`nwbinspector`](https://nwbinspector.readthedocs.io/) (CLI)
- [`dandi`](https://www.dandiarchive.org/handbook/13_upload/) (CLI)

Install with pip:

```bash
pip install pynwb hdmf pynapple pandas scipy pytz openpyxl nwbinspector dandi
```

## Configuration

Edit the constants near the top of `PeyracheToNWB.py`:

```python
ADD_RAW = False              # True to include the raw 20 kHz .dat (large)
OVERWRITE_EXISTING = True    # False to skip sessions that already have a .nwb
DATAPATH = Path('/Volumes/Extreme SSD/Dataset_Main_PoSub')
DANDIPATH = DATAPATH / 'NWB' / '000939'
METADATA_FILENAME = 'Dataset_metadata.xlsx'
```

## Running

```bash
python PeyracheToNWB.py
```

The script iterates every row in `Dataset_metadata.xlsx`, writes `<session>.nwb` next to the source data, validates it with `nwbinspector --config dandi`, and then runs `dandi organize` into `DANDIPATH`.

## Citing

If you use this dataset, please cite the paper:

> Duszkiewicz AJ, *et al.* *Nature Neuroscience* (2024). [doi:10.1038/s41593-024-01588-5](https://doi.org/10.1038/s41593-024-01588-5)

## License

See [LICENSE](LICENSE).
