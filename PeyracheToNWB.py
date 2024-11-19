
import os
import warnings
import subprocess
import math
from datetime import datetime
import pytz
import pandas as pd
import pynapple as nap
from matplotlib.pyplot import *
import scipy.io as spio
from pathlib import Path
from pynwb import NWBFile
from pynwb import NWBHDF5IO
from pynwb.file import Subject
from pynwb.epoch import TimeIntervals
from pynwb.behavior import Position
from pynwb.behavior import CompassDirection
from pynwb.ecephys import TimeSeries
from pynwb.ecephys import ElectricalSeries
from pynwb.ecephys import LFP
from pynwb.ogen import OptogeneticSeries
from pynwb.ogen import OptogeneticStimulusSite
from hdmf.backends.hdf5.h5_utils import H5DataIO
from pynwb.behavior import SpatialSeries


# TO DO:

# Delete NWB file if there is one already

# Before creating NWB, check:
# mouse strain, virus type, epochs, opto device (where was it implanted), ephys device (angle of implant)

# tutorial: https://www.youtube.com/watch?v=rlywed3ar-s&ab_channel=NeurodataWithoutBorders
# doc: https://nwb-overview.readthedocs.io/


addRaw = False
removeNWB = True  # False to skip NWBs already present, True to overwrite

datapath = Path('/Volumes/Extreme SSD/Dataset_Main_PoSub')
dandipath = datapath / 'NWB' / '000939'
os.chdir(datapath)

# read metadata to get the number of folders
print('Reading metadata from the xml file... ')
metaname = 'Dataset_metadata.xlsx'
metadata_full = pd.read_excel(datapath / metaname)
totFolders = len(metadata_full)
folder_ids = range(totFolders)

def main():
    for nfolder in folder_ids:

        foldername = metadata_full.iloc[nfolder]
        foldername = foldername['Recording']
        print('Converting folder ' + foldername + '')

        make_nwb(nfolder)
        organize_nwb(nfolder)

def organize_nwb(rec_number):

    # first run nwbinspector
    os.chdir(datapath)
    metadata = metadata_full.iloc[rec_number]
    foldername = metadata['Recording']
    path = datapath / foldername
    filename = foldername + '.nwb'

    # first run nwbinspector
    print('Inspecting NWB file...')
    command = ['nwbinspector', foldername,  '--config', 'dandi']
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    if 'Found 0 issues' not in result.stdout:
        print('INSPECTION FAILED :(')
        return
    
    # now organize the NWB file into the
    print('Moving to dandiset...')
    os.chdir(dandipath)
    command = ['dandi', 'organize', str(path)]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print(result.stdout)
    except Exception as e:
        print(f"An error occurred: {e}")
        return
        
    print('Done!')


def make_nwb(rec_number):

    # pick correct metadata
    metadata = metadata_full.iloc[rec_number]
    foldername = metadata['Recording']
    path = datapath / foldername / 'Data'

    # remove or skip NWB file if one is present in folder
    nwbpath = datapath / foldername / str(foldername + '.nwb')
    if os.path.exists(str(nwbpath)):
        if removeNWB:
            os.remove(nwbpath)
            print('Old NWB file removed')
        else:
            return

    # get info from folder name
    fileInfo = foldername.split('-')
    start_time = datetime(int(('20' + fileInfo[1][0:2])), int((fileInfo[1][2:4])), int((fileInfo[1][4:6])))
    start_time = pytz.timezone('America/New_York').localize(start_time)

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

   # add subject
    nwbfile.subject = Subject(age='P12W/',
                              description=metadata['Mouse_line'],
                              species='Mus musculus',
                              subject_id=fileInfo[0],
                              sex='M')


    # Load tracking, epochs and spikes
    print('Loading variables from .mat files...')
    pos, ang, epochs, spikes, shank_id, waveforms, maxIx, tr2pk = import_session(path)

    # loading channel order and good channels
    chanmap_file = datapath / foldername / 'ChannelMap.mat'
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
    nwbfile.add_epoch(start_time=epochs['Start'][0], stop_time=epochs['End'][0], tags=metadata['Epoch_1'])
    nwbfile.add_epoch(start_time=epochs['Start'][1], stop_time=epochs['End'][1], tags=metadata['Epoch_2'])
    if not isinstance(metadata['Epoch_3'], float):
        nwbfile.add_epoch(start_time=epochs['Start'][2], stop_time=epochs['End'][2], tags=metadata['Epoch_3'])
    if not isinstance(metadata['Epoch_3'], float):
        nwbfile.add_epoch(start_time=epochs['Start'][3], stop_time=epochs['End'][3], tags=metadata['Epoch_4'])

    sleep_file = datapath/foldername/'Sleep'/(foldername + '.SleepState.states.mat')
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
    acc_file = datapath / foldername / (foldername + '_auxiliary.dat')
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
    # set up shank(s)
    nwbfile.add_electrode_column(name='label', description='label of electrode')
    nwbfile.add_electrode_column(name='is_faulty', description='Boolean column to indicate faulty electrodes')
    step = 12.5  # difference in spacing between electrodes

    n_shanks = 1
    n_channels = 64

    electrode_counter = 0
    device = nwbfile.create_device(
        name='Cambridge Neurotech H5 probe',
        description=metadata['Probe_description']
    )

    for ishank in range(n_shanks):
        #create an electrode group for this shank
        electrode_group = nwbfile.create_electrode_group(
            name='shank{}'.format(ishank),
            description='electrode group for shank{}'.format(ishank),
            device=device,
            location='Postsubiculum (left hemisphere)',
        )

        for ielec in range(n_channels):
            elec_depth = step * (n_channels - ielec-1)
            nwbfile.add_electrode(
                x=0., y=elec_depth, z=0.,  # add electrode position
                location='Postsubiculum (left hemisphere)',
                filtering='none',
                is_faulty=isFaulty[electrode_counter],
                group=electrode_group,
                label='shank{}elec{}'.format(ishank,ielec)
            )
            electrode_counter += 1

    # define table region RAW DAT FILE and LFP will refer to (all electrodes)
    all_table_region = nwbfile.create_electrode_table_region(
        region=list(range(electrode_counter)),
        description='all electrodes',
    )

    # RAW DAT FILE
    if addRaw:
        print('Adding raw dat file (may take a wee while)...')
        path_raw = datapath / foldername / (foldername + '.dat')
        raw_data = nap.load_eeg(filepath=path_raw, channel=None, n_channels=64, frequency=20000.0, precision='int16',
                                bytes_size=2)
        raw_data = raw_data[:, chanOrder]  # sort according to channel order

        raw_electrical_series = ElectricalSeries(
            name="ElectricalSeries",
            data=H5DataIO(raw_data, compression=True),  # use this function to compress
            electrodes=all_table_region,
            starting_time=0.0,  # timestamp of the first sample in seconds relative to the session start time
            rate=20000.0,  # in Hz
        )

        nwbfile.add_acquisition(raw_electrical_series)

    # LFP
    print('Adding lfp...')

    path_lfp = datapath / foldername / (foldername + '.lfp')
    lfp_data = nap.load_eeg(filepath=path_lfp, channel=None, n_channels=64, frequency=1250.0, precision='int16',
                       bytes_size=2)
    lfp_data = lfp_data[:, chanOrder]  # sort according to channel order


    # create ElectricalSeries
    lfp_elec_series = ElectricalSeries(
        name='LFP',
        data=H5DataIO(lfp_data, compression=True),  # use this function to compress
        description='Local field potential (low-pass filtered at 625 Hz)',
        electrodes=all_table_region,
        rate=1250.
    )

    # store ElectricalSeries in an LFP container
    warnings.filterwarnings("ignore", message=".*DynamicTableRegion.*")  # this is to supress a warning here that doesn't seem cause any issues
    lfp = LFP(electrical_series=lfp_elec_series)
    warnings.resetwarnings()

    ecephys_module = nwbfile.create_processing_module(name='ecephys',
                                                      description='Processed electrophysiological signals'
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

    for ncell in range(len(spikes)):
        nwbfile.add_unit(spike_times=spikes[ncell].times(),
                         electrode_index=maxIx[ncell]-1,
                         waveform_mean=waveforms[ncell].T,
                         trough_to_peak=tr2pk[ncell],
                         is_excitatory=isEX[ncell],
                         is_fast_spiking=isFS[ncell],
                         is_head_direction=isHD[ncell],
                         electrode_group=nwbfile.electrode_groups['shank0'])

    # EMG
    print('Adding emg...')
    emg_file = datapath / foldername / 'Sleep' / (foldername + '.EMGFromLFP.LFP.mat')
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
    opto_file = datapath / foldername / 'Opto_digitalin.dat'
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

    with NWBHDF5IO(datapath/foldername/(foldername + '.nwb'), 'w') as io:
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

    # Next lines load the spike data from the .mat file
    spikedata = spio.loadmat(spike_file, simplify_cells=True)
    total_cells = np.arange(0, len(spikedata['S']['C']))  # To find the total number of cells in the recording
    spikes = dict()  # For pynapple, this will be turned into a TsGroup of all cells' timestamps
    cell_df = pd.DataFrame(columns=['timestamps'])  # Dataframe for cells and their spike timestamps
    # Loop to assign cell timestamps into the dataframe and the cell_ts dictionary
    for cell in total_cells:  # give spikes
        timestamps = spikedata['S']['C'][cell]['tsd']['t']
        cell_df.loc[cell, 'timestamps'] = timestamps
        temp = {cell: nap.Ts(timestamps)}
        spikes.update(temp)

   # get shank and channel ID for cells
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