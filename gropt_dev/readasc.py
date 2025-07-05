"""
Read ASC files from Siemens scanners.
Taken directly from pypulseq: https://github.com/imr-framework/pypulseq/
"""

import re
from types import SimpleNamespace
from typing import List, Tuple

import numpy as np


def readasc(filename: str) -> Tuple[dict, dict]:
    """
    Read Siemens ASC ascii-formatted textfile and returns a dictionary structure.
    
    E.g. a[0].b[2][3].c = "string"
    parses into:
      asc['a'][0]['b'][2][3]['c'] = "string"

    Parameters
    ----------
    filename : str
        Filename of the ASC file.

    Returns
    -------
    asc : dict
        Dictionary of ASC part of file.
    extra : dict
        Dictionary of other fields after "ASCCONV END"
    """
    asc, extra = {}, {}

    # Read asc file and convert it into a dictionary structure
    with open(filename, 'r') as fp:
        end_of_asc = False

        for next_line in fp:
            next_line = next_line.strip()

            if next_line == '### ASCCONV END ###':  # find end of mrProt in the asc file
                end_of_asc = True

            if next_line == '' or next_line[0] == '#':
                continue

            # regex wizardry: Matches lines like 'a[0].b[2][3].c = "string" # comment'
            # Note this assumes correct formatting, e.g. does not check whether
            # brackets match.
            match = re.match(
                r'^\s*([a-zA-Z0-9\[\]\._]+)\s*\=\s*(("[^"]*"|\'[^\']\')|(\d+)|([0-9\.e\-]+))\s*((#|\/\/)(.*))?$',
                next_line,
            )

            if match:
                field_name = match[1]

                # Keep track of where to put the value: base[assign_to] = value
                if end_of_asc:
                    base = extra
                else:
                    base = asc

                assign_to = None

                # Iterate over every segment of the field name
                parts = field_name.split('.')
                for p in parts:
                    # Update base so final assignment is like: base[assign_to][p] = value
                    if assign_to is not None and assign_to not in base:
                        base[assign_to] = {}
                    if assign_to is not None:
                        base = base[assign_to]

                    # Iterate over brackets
                    start = p.find('[')
                    if start != -1:
                        name = p[:start]
                        assign_to = name

                        while start != -1:
                            stop = p.find(']', start)
                            index = int(p[start + 1 : stop])

                            # Update base so final assignment is like: base[assign_to][p][index] = value
                            if assign_to not in base:
                                base[assign_to] = {}
                            base = base[assign_to]
                            assign_to = index

                            start = p.find('[', stop)
                    else:
                        assign_to = p

                # Depending on which regex section matched we can infer the value type
                if match[3]:
                    base[assign_to] = match[3][1:-1]
                elif match[4]:
                    base[assign_to] = int(match[4])
                elif match[5]:
                    base[assign_to] = float(match[5])
                else:
                    raise RuntimeError('This should not be reached')
            elif next_line.find('=') != -1:
                raise RuntimeError(f'Bug: ASC line with an assignment was not parsed correctly: {next_line}')

    return asc, extra



def asc_to_acoustic_resonances(asc: dict) -> List[dict]:
    """
    Convert ASC dictionary from readasc to list of acoustic resonances.

    Parameters
    ----------
    asc : dict
        ASC dictionary, see readasc

    Returns
    -------
    List[dict]
        List of acoustic resonances (specified by frequency and bandwidth fields).
    """
    if 'aflGCAcousticResonanceFrequency' in asc:
        freqs = asc['aflGCAcousticResonanceFrequency']
        bw = asc['aflGCAcousticResonanceBandwidth']
    else:
        freqs = asc['asGPAParameters'][0]['sGCParameters']['aflAcousticResonanceFrequency']
        bw = asc['asGPAParameters'][0]['sGCParameters']['aflAcousticResonanceBandwidth']

    return [{'frequency': f, 'bandwidth': b} for f, b in zip(freqs.values(), bw.values()) if f != 0]


def asc_to_hw(asc: dict, cardiac_model: bool = False) -> SimpleNamespace:
    """
    Convert ASC dictionary from readasc to SAFE hardware description.

    Parameters
    ----------
    asc : dict
        ASC dictionary, see readasc
    cardiac_model : bool
        Whether or not to read the cardiac stimulation model instead of the
        default PNS model (returns None if not available)

    Returns
    -------
    SimpleNamespace
        SAFE hardware description
    """
    hw = SimpleNamespace()

    if 'asCOMP' in asc and 'tName' in asc['asCOMP']:
        hw.name = asc['asCOMP']['tName']
    else:
        hw.name = 'unknown'

    if 'GradPatSup' in asc:
        asc_pns = asc['GradPatSup']['Phys']['PNS']
    else:
        asc_pns = asc

    if cardiac_model:
        if 'GradPatSup' in asc and 'CarNS' in asc['GradPatSup']['Phys']:
            asc_pns = asc['GradPatSup']['Phys']['CarNS']
        else:
            return None

    hw.x = SimpleNamespace()
    hw.x.tau1 = asc_pns['flGSWDTauX'][0]  # ms
    hw.x.tau2 = asc_pns['flGSWDTauX'][1]  # ms
    hw.x.tau3 = asc_pns['flGSWDTauX'][2]  # ms
    hw.x.a1 = asc_pns['flGSWDAX'][0]
    hw.x.a2 = asc_pns['flGSWDAX'][1]
    hw.x.a3 = asc_pns['flGSWDAX'][2]
    hw.x.stim_limit = asc_pns['flGSWDStimulationLimitX']  # T/m/s
    hw.x.stim_thresh = asc_pns['flGSWDStimulationThresholdX']  # T/m/s

    hw.y = SimpleNamespace()
    hw.y.tau1 = asc_pns['flGSWDTauY'][0]  # ms
    hw.y.tau2 = asc_pns['flGSWDTauY'][1]  # ms
    hw.y.tau3 = asc_pns['flGSWDTauY'][2]  # ms
    hw.y.a1 = asc_pns['flGSWDAY'][0]
    hw.y.a2 = asc_pns['flGSWDAY'][1]
    hw.y.a3 = asc_pns['flGSWDAY'][2]
    hw.y.stim_limit = asc_pns['flGSWDStimulationLimitY']  # T/m/s
    hw.y.stim_thresh = asc_pns['flGSWDStimulationThresholdY']  # T/m/s

    hw.z = SimpleNamespace()
    hw.z.tau1 = asc_pns['flGSWDTauZ'][0]  # ms
    hw.z.tau2 = asc_pns['flGSWDTauZ'][1]  # ms
    hw.z.tau3 = asc_pns['flGSWDTauZ'][2]  # ms
    hw.z.a1 = asc_pns['flGSWDAZ'][0]
    hw.z.a2 = asc_pns['flGSWDAZ'][1]
    hw.z.a3 = asc_pns['flGSWDAZ'][2]
    hw.z.stim_limit = asc_pns['flGSWDStimulationLimitZ']  # T/m/s
    hw.z.stim_thresh = asc_pns['flGSWDStimulationThresholdZ']  # T/m/s

    if 'asGPAParameters' in asc:
        hw.x.g_scale = asc['asGPAParameters'][0]['sGCParameters']['flGScaleFactorX']
        hw.y.g_scale = asc['asGPAParameters'][0]['sGCParameters']['flGScaleFactorY']
        hw.z.g_scale = asc['asGPAParameters'][0]['sGCParameters']['flGScaleFactorZ']
    else:
        print('Warning: Gradient scale factors not in ASC file: assuming 1/pi')
        hw.x.g_scale = 1 / np.pi
        hw.y.g_scale = 1 / np.pi
        hw.z.g_scale = 1 / np.pi

    return hw

def hw_namespace_to_dict(hw: SimpleNamespace) -> dict:
    """
    Convert a SimpleNamespace object representing hardware parameters to a dictionary.

    Also converts tau parameters from milliseconds to seconds.
    
    Parameters
    ----------
    hw : SimpleNamespace
        Hardware parameters from pyPulseq asc_to_hw function.

    Returns
    -------
    dict
        Dictionary representation of hardware parameters in seconds.
    """
    out = {
        'tau1': np.array([hw.x.tau1, hw.y.tau1, hw.z.tau1])/1000.0,  # Convert ms to seconds
        'tau2': np.array([hw.x.tau2, hw.y.tau2, hw.z.tau2])/1000.0,  # Convert ms to seconds
        'tau3': np.array([hw.x.tau3, hw.y.tau3, hw.z.tau3])/1000.0,  # Convert ms to seconds
        'a1': np.array([hw.x.a1, hw.y.a1, hw.z.a1]),
        'a2': np.array([hw.x.a2, hw.y.a2, hw.z.a2]),
        'a3': np.array([hw.x.a3, hw.y.a3, hw.z.a3]),
        'stim_limit': np.array([hw.x.stim_limit, hw.y.stim_limit, hw.z.stim_limit]),
        'g_scale': np.array([hw.x.g_scale, hw.y.g_scale, hw.z.g_scale]),
    }
        
    return out

def asc_to_safe(asc_file: str) -> tuple[dict, dict]:
    """
    Convert an ASC file to SAFE format.

    Parameters
    ----------
    asc_file : str
        Path to the ASC file.

    Returns
    -------
    tuple[dict, dict]
        SAFE hardware parameters. The first dictionary contains the parameters for the
        PNS model, and the second dictionary contains the parameters for the cardiac model.
        If the cardiac model is not available, the second dictionary will be None.
    """
    asc, extra = readasc(asc_file)

    try:
        hw = asc_to_hw(asc, cardiac_model=False)
        hw_cardiac = asc_to_hw(asc, cardiac_model=True)
    except KeyError as e:
        print(f'ERROR: finding SAFE parameters in asc file. KeyError: {e}')
        hw = None
        hw_cardiac = None
        
    if hw is not None:
        hw = hw_namespace_to_dict(hw)
    if hw_cardiac is not None:
        hw_cardiac = hw_namespace_to_dict(hw_cardiac)
        
    return hw, hw_cardiac
