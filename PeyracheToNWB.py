
import os
import math
import warnings
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytz
import pynapple as nap
import scipy.io as spio
from hdmf.backends.hdf5.h5_utils import H5DataIO
from pynwb import NWBFile, NWBHDF5IO
from pynwb.behavior import CompassDirection, Position, SpatialSeries
from pynwb.ecephys import ElectricalSeries, LFP, TimeSeries
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from pynwb.ogen import OptogeneticSeries, OptogeneticStimulusSite


# TO DO:

# Delete NWB file if there is one already

# Before creating NWB, check:
# mouse strain, virus type, epochs, opto device (where was it implanted), ephys device (angle of implant)

# tutorial: https://www.youtube.com/watch?v=rlywed3ar-s&ab_channel=NeurodataWithoutBorders
# doc: https://nwb-overview.readthedocs.io/


ADD_RAW = False
OVERWRITE_EXISTING = True  # False to skip NWBs already present, True to overwrite

DATAPATH = Path('/Volumes/Extreme SSD/Dataset_Main_PoSub')
DANDIPATH = DATAPATH / 'NWB' / '000939'
METADATA_FILENAME = 'Dataset_metadata.xlsx'


def _get_meta(metadata, key, default=None):
    """Return metadata[key] if present and not NaN, otherwise default.

    Lets the metadata spreadsheet grow new columns without breaking older
    sessions that don't have them.
    """
    if key not in metadata.index:
        return default
    val = metadata[key]
    if isinstance(val, float) and math.isnan(val):
        return default
    return val


def main():
    print('Reading metadata from the Excel file...')
    metadata_full = pd.read_excel(DATAPATH / METADATA_FILENAME)

    for nfolder in range(len(metadata_full)):
        foldername = metadata_full.iloc[nfolder]['Recording']
        print(f'Converting folder {foldername}')

        make_nwb(metadata_full, nfolder)
        organize_nwb(metadata_full, nfolder)

def organize_nwb(metadata_full, rec_number):
    metadata = metadata_full.iloc[rec_number]
    foldername = metadata['Recording']
    session_path = DATAPATH / foldername

    print('Inspecting NWB file...')
    result = subprocess.run(
        ['nwbinspector', foldername, '--config', 'dandi'],
        capture_output=True, text=True, cwd=DATAPATH,
    )
    if result.returncode != 0 or 'Found 0 issues' not in result.stdout:
        print('INSPECTION FAILED :(')
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
        return

    print('Moving to dandiset...')
    result = subprocess.run(
        ['dandi', 'organize', str(session_path)],
        capture_output=True, text=True, cwd=DANDIPATH,
    )
    print(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr)
        return

    print('Done!')


def make_nwb(metadata_full, rec_number):

    # pick correct metadata
    metadata = metadata_full.iloc[rec_number]
    foldername = metadata['Recording']
    path = DATAPATH / foldername / 'Data'

    # remove or skip NWB file if one is present in folder
    nwbpath = DATAPATH / foldername / f'{foldername}.nwb'
    if nwbpath.exists():
        if OVERWRITE_EXISTING:
            nwbpath.unlink()
            print('Old NWB file removed')
        else:
            return

    # get info from folder name (expected format: <subject_id>-YYMMDD-...)
    fileInfo = foldername.split('-')
    start_time = datetime.strptime(fileInfo[1], '%y%m%d')
    start_time = pytz.timezone('America/Toronto').localize(start_time)

    # create an nwb file
    print('Creating NWB file and adding metadata...')
    nwbfile = NWBFile(
        session_description='Open field and sleep recording',
        experiment_description='high-density recordings in mouse postsubiculum during exploration and sleep',
        identifier=fileInfo[0],
        session_start_time=start_time,
        session_id=fileInfo[1],
        experimenter='Duszkiewicz, Adrian J.',
        lab='Peyrache Lab',
        institution='McGill University',
        virus=str(metadata['Virus']),
        related_publications='doi: 10.1038/s41593-024-01588-5',
        keywords=['head-direction', 'postsubiculum', 'extracellular', 'freely-moving', 'electrophysiology']

    )

    # add subject (age/sex come from metadata if present)
    nwbfile.subject = Subject(
        age=_get_meta(metadata, 'Age', default='P12W'),
        description=metadata['Mouse_line'],
        species='Mus musculus',
        subject_id=fileInfo[0],
        sex=_get_meta(metadata, 'Sex', default='M'),
    )


    # Load tracking, epochs and spikes
    print('Loading variables from .mat files...')
    pos, ang, epochs, spikes, shank_id, waveforms, maxIx, tr2pk = import_session(path)

    # loading channel order and good channels
    chanmap_file = DATAPATH / foldername / 'ChannelMap.mat'
    chanmap = spio.loadmat(chanmap_file, simplify_cells=True)
    chanOrder = chanmap['chanOrder']
    goodChans = chanmap['goodChans'].astype(bool)
    # get good channels
    isFaulty = ~goodChans[chanOrder]  # sort according to channel order


    #### BEHAVIOUR #####

    # create behaviour module
    behavior_module = nwbfile.create_processing_module(
        name='behavior',
        description='Tracking data acquired by Optitrack Motive 2.0'
    )

    # EPOCHS
    print('Adding epochs and behavioural variables...')
    for epoch_idx, epoch_key in enumerate(['Epoch_1', 'Epoch_2', 'Epoch_3', 'Epoch_4']):
        tag = _get_meta(metadata, epoch_key)
        if tag is None or epoch_idx >= len(epochs):
            continue
        nwbfile.add_epoch(
            start_time=epochs['Start'][epoch_idx],
            stop_time=epochs['End'][epoch_idx],
            tags=tag,
        )

    sleep_file = DATAPATH/foldername/'Sleep'/(foldername + '.SleepState.states.mat')
    sleepEpochs = spio.loadmat(sleep_file, simplify_cells=True)
    epWake = np.float32(sleepEpochs['SleepState']['ints']['WAKEstate'])
    epNREM = np.float32(sleepEpochs['SleepState']['ints']['NREMstate'])
    epREM = np.float32(sleepEpochs['SleepState']['ints']['REMstate'])


    if epREM.size > 0:
        rem = TimeIntervals(name='rem')
        if epREM.ndim == 1:  # in case there is only one interval
            rem.add_row(start_time=epREM[0], stop_time=epREM[1], tags=str(0))
        elif epREM.ndim == 2:
            for nrow in range(len(epREM)):
                rem.add_row(start_time=epREM[nrow, 0], stop_time=epREM[nrow, 1], tags=str(nrow))
        nwbfile.add_time_intervals(rem)

    if epNREM.size > 0:
        nrem = TimeIntervals(name='nrem')
        if epNREM.ndim == 1: # in case there is only one interval
            nrem.add_row(start_time=epNREM[0], stop_time=epNREM[1], tags=str(0))
        elif epNREM.ndim == 2:
            for nrow in range(len(epNREM)):
                nrem.add_row(start_time=epNREM[nrow, 0], stop_time=epNREM[nrow, 1], tags=str(nrow))
        nwbfile.add_time_intervals(nrem)

    # TRACKING
    spatial_series_obj = SpatialSeries(
        name='position',
        description='(x,y) position in the open field',
        data=pos.values,
        timestamps=pos.index.to_numpy(),
        reference_frame='Cue is on the top wall (highest value of y)',
        unit='centimeters'
    )
    position_obj = Position(spatial_series=spatial_series_obj)

    spatial_series_obj = SpatialSeries(
        name='head-direction',
        description='Horizontal angle of the head (yaw)',
        data=ang.values,
        timestamps=ang.index.to_numpy(),
        reference_frame='Clockwise, 0 rad = top (high y values)',
        unit='radians'
    )
    direction_obj = CompassDirection(spatial_series=spatial_series_obj)

    # update behaviour module
    behavior_module.add(position_obj)
    behavior_module.add(direction_obj)

    # ACCELEROMETER
    acc_file = DATAPATH / foldername / (foldername + '_auxiliary.dat')
    if os.path.exists(acc_file):
        print('Adding accelerometer data...')
        acc_data = nap.load_eeg(filepath=acc_file, channel=None, n_channels=3, frequency=20000, precision='int16',
                                bytes_size=2)

        acc = TimeSeries(
            name="accelerometer",
            description="Accelerometer data from the Intan RHD headstage",
            data=H5DataIO(acc_data, compression=True),  # use this function to compress
            #data=acc_data,
            unit="a.u.",
            rate=20000.,
        )

        behavior_module.add(acc)
    else:
        print('Accelerometer data not found.')

    ### EPHYS ###

    print('Adding electrodes...')
    nwbfile.add_electrode_column(name='label', description='label of electrode')
    nwbfile.add_electrode_column(name='is_faulty', description='Boolean column to indicate faulty electrodes')
    step = 12.5  # difference in spacing between electrodes

    # Derive geometry/sample rates from the channel map and metadata so the
    # script works on recordings with different probe / acquisition configs.
    n_channels = int(len(chanOrder))
    n_shanks = int(_get_meta(metadata, 'N_shanks', default=1))
    raw_rate = float(_get_meta(metadata, 'Raw_sample_rate', default=20000.0))
    lfp_rate = float(_get_meta(metadata, 'LFP_sample_rate', default=1250.0))
    if n_channels % n_shanks != 0:
        raise ValueError(
            f'n_channels ({n_channels}) is not divisible by n_shanks ({n_shanks})'
        )
    channels_per_shank = n_channels // n_shanks

    device = nwbfile.create_device(
        name='Cambridge Neurotech H5 probe',
        description=metadata['Probe_description']
    )

    electrode_counter = 0
    for ishank in range(n_shanks):
        electrode_group = nwbfile.create_electrode_group(
            name=f'shank{ishank}',
            description=f'electrode group for shank{ishank}',
            device=device,
            location='Postsubiculum (left hemisphere)',
        )

        for ielec in range(channels_per_shank):
            elec_depth = step * (channels_per_shank - ielec - 1)
            nwbfile.add_electrode(
                x=0., y=elec_depth, z=0.,
                location='Postsubiculum (left hemisphere)',
                filtering='none',
                is_faulty=bool(isFaulty[electrode_counter]),
                group=electrode_group,
                label=f'shank{ishank}elec{ielec}',
            )
            electrode_counter += 1

    # define table region RAW DAT FILE and LFP will refer to (all electrodes)
    all_table_region = nwbfile.create_electrode_table_region(
        region=list(range(electrode_counter)),
        description='all electrodes',
    )

    # RAW DAT FILE
    if ADD_RAW:
        print('Adding raw dat file (may take a wee while)...')
        path_raw = DATAPATH / foldername / (foldername + '.dat')
        raw_data = nap.load_eeg(filepath=path_raw, channel=None, n_channels=n_channels,
                                frequency=raw_rate, precision='int16', bytes_size=2)
        raw_data = raw_data[:, chanOrder]  # sort according to channel order

        raw_electrical_series = ElectricalSeries(
            name="ElectricalSeries",
            data=H5DataIO(raw_data, compression=True),
            electrodes=all_table_region,
            starting_time=0.0,
            rate=raw_rate,
        )

        nwbfile.add_acquisition(raw_electrical_series)

    # LFP
    print('Adding lfp...')

    path_lfp = DATAPATH / foldername / (foldername + '.lfp')
    lfp_data = nap.load_eeg(filepath=path_lfp, channel=None, n_channels=n_channels,
                            frequency=lfp_rate, precision='int16', bytes_size=2)
    lfp_data = lfp_data[:, chanOrder]  # sort according to channel order

    lfp_elec_series = ElectricalSeries(
        name='LFP',
        data=H5DataIO(lfp_data, compression=True),
        description='Local field potential (low-pass filtered at 625 Hz)',
        electrodes=all_table_region,
        rate=lfp_rate,
    )

    # Scope the warning suppression so we don't clobber other filters.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*DynamicTableRegion.*")
        lfp = LFP(electrical_series=lfp_elec_series)

    ecephys_module = nwbfile.create_processing_module(
        name='ecephys',
        description='Processed electrophysiological signals',
    )
    ecephys_module.add(lfp)

    # UNITS
    print('adding spikes...')
    nwbfile.add_unit_column(name="electrode_index", description="electrode with the highest waveform amplitude")
    nwbfile.add_unit_column(name="trough_to_peak", description="Trough-to-peak duration of waveform (ms)")
    nwbfile.add_unit_column(name="is_excitatory", description="Putative excitatory cells")
    nwbfile.add_unit_column(name="is_fast_spiking", description="Putative fast-spiking interneurons")
    nwbfile.add_unit_column(name="is_head_direction", description="Head-direction tuned units")

    #load cell types
    celltypes_file = path / 'CellTypes.mat'
    celltypes = spio.loadmat(celltypes_file, simplify_cells=True)
    isHD = celltypes['hd']
    isEX = celltypes['ex']
    isFS = celltypes['fs']

    # shank_id is 1-indexed in MATLAB; maxIx is also 1-indexed and waveforms
    # are stored as (channels, samples) and need transposing to NWB's
    # (samples, channels) convention.
    shank_id_arr = np.atleast_1d(shank_id).astype(int)
    for ncell in range(len(spikes)):
        unit_shank = int(shank_id_arr[ncell]) - 1 if shank_id_arr.size else 0
        group_name = f'shank{unit_shank}'
        if group_name not in nwbfile.electrode_groups:
            group_name = 'shank0'  # fall back if cluster shank > probe shanks
        nwbfile.add_unit(
            spike_times=spikes[ncell].times(),
            electrode_index=int(maxIx[ncell]) - 1,
            waveform_mean=waveforms[ncell].T,
            trough_to_peak=tr2pk[ncell],
            is_excitatory=isEX[ncell],
            is_fast_spiking=isFS[ncell],
            is_head_direction=isHD[ncell],
            electrode_group=nwbfile.electrode_groups[group_name],
        )

    # EMG
    print('Adding emg...')
    emg_file = DATAPATH / foldername / 'Sleep' / (foldername + '.EMGFromLFP.LFP.mat')
    emg = spio.loadmat(emg_file, simplify_cells=True)
    emg = emg['EMGFromLFP']['data']

    emg = TimeSeries(
        name="pseudoEMG",
        description="Pseudo EMG from correlated high-frequency LFP",
        data=emg,
        unit="a.u.",
        starting_time=0.,
        rate=2.0,
    )

    ecephys_module.add(emg)

    # OPTO
    opto_file = DATAPATH / foldername / 'Opto_digitalin.dat'
    if os.path.exists(opto_file):
        print('Adding optogenetic series...')
        opto_device = nwbfile.create_device(
            name='Mono fiberoptic cannula',
            description=metadata['Opto_device'],
            manufacturer='Doric Lenses'
        )
        opto_site = OptogeneticStimulusSite(
            name=metadata['Opto_site_name'],
            device=opto_device,
            description=metadata['Opto_site_description'],
            excitation_lambda=metadata['Wavelength'],
            location=metadata['Opto_site_location']
        )
        nwbfile.add_ogen_site(opto_site)
        ogen_module = nwbfile.create_processing_module(name='ogen',
                                                       description='Optogenetic stimulation data')

        opto_data = nap.load_eeg(filepath=opto_file, channel=None, n_channels=1, frequency=20000.0, precision='int16',
                                bytes_size=2)

        opto = OptogeneticSeries(
            name='optogenetic_stim',
            description='optogenetic stimulation',
            site=opto_site,
            conversion=1.,
            rate=20000.0,
            data=H5DataIO(opto_data.ravel(), compression=True),  # use this function to compress (ravel to convert to 1D)
            starting_time=0.,
        )

        ogen_module.add(opto)
    else:
        print('No opto file found.')


    # save NWB file
    print('Saving NWB file...')

    with NWBHDF5IO(DATAPATH/foldername/(foldername + '.nwb'), 'w') as io:
        io.write(nwbfile)

    print('Done!')



def import_session(path):
    # file names for spikes, angle, and epoch files
    spike_file = path / 'SpikeData.mat'
    angle_file = path / 'Angle.mat'  # has HD angle
    epoch_file = path / 'Epoch_TS.csv'
    pos_file = path / 'Position.mat'
    waveform_file = path / 'Waveforms.mat'
    wfeatures_file = path / 'WaveformFeatures.mat'

    # load epochs from mat file
    epochs = pd.read_csv(epoch_file, header=None, names=['Start','End'])

    # Load angle from mat file
    angle_data = spio.loadmat(angle_file, simplify_cells=True)
    ang = angle_data['ang']['data']
    ang = ang % (2 * np.pi)  # mod 2 pi
    ang = pd.Series(ang)
    ang.index = angle_data['ang']['t']

    # load position from mat file
    pos_data = spio.loadmat(pos_file, simplify_cells=True)
    pos = pd.DataFrame()
    pos['X'] = pos_data['pos']['data'][:, 0]
    pos['Y'] = pos_data['pos']['data'][:, 1]
    pos.index = pos_data['pos']['t']

    # Load spike timestamps per cell
    spikedata = spio.loadmat(spike_file, simplify_cells=True)
    spikes = {
        cell: nap.Ts(spikedata['S']['C'][cell]['tsd']['t'])
        for cell in range(len(spikedata['S']['C']))
    }
    shank_id = spikedata['shank']

    # get waveforms and waveform features
    waveforms = spio.loadmat(waveform_file, simplify_cells=True)
    waveforms = waveforms['meanWaveforms']
    wfeatures = spio.loadmat(wfeatures_file, simplify_cells=True)
    maxIx = wfeatures['maxIx']
    tr2pk = wfeatures['tr2pk']


    return pos, ang, epochs, spikes, shank_id, waveforms, maxIx, tr2pk


if __name__ == '__main__':
    main()