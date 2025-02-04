import os
from os import listdir, path
from ev_model import elements, environment
from lxml import etree
import pandas as pd
import sys

def save_0xml_file(evs, n, dt, degrees_of_freedom, base_filename,
        target_source_points_per_cell, base_dir='resources/iterations', version=0,
        boundaries=False, ev_ev_collisions=False, cells=None, new_evs=False,
        new_evs_interval_seconds=30.0, new_evs_threshold=0.9,
        seconds_in_initial_state=2.0, min_ev_radius=40, max_ev_radius=120,
        ev_viability=3600):

    experiment_directory = f"{base_filename}_{n/1000:.1f}k_{'intCol' if ev_ev_collisions else 'noCol'}{'_newEvs' if new_evs else ''}"
    experiment_name = f"{n/1000:.1f}k_{'intCol' if ev_ev_collisions else 'noCol'}{'_newEVs' if new_evs else ''}"

    f0 = f"{base_dir}/v{version}/{experiment_directory}/0_{experiment_name}.xml"

    directory = os.path.dirname(f0)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(f0, 'w') as of:
        of.write(f"""<states>
    <itno>0</itno>
    <environment>
        <ev_collisions>{1 if ev_ev_collisions else 0}</ev_collisions>
        <boundaries_enabled>{1 if boundaries and cells else 0}</boundaries_enabled>
        <introduce_new_evs>{1 if new_evs else 0}</introduce_new_evs>
        <seconds_before_introducing_new_evs>{new_evs_interval_seconds if new_evs else 100000 :.1f}</seconds_before_introducing_new_evs>
        <new_evs_threshold>{new_evs_threshold if new_evs else 1 :.2f}</new_evs_threshold>
        <dt>{dt}</dt>
        <dof>{degrees_of_freedom}</dof>
        <min_ev_radius>{min_ev_radius}</min_ev_radius>
        <max_ev_radius>{max_ev_radius}</max_ev_radius>
        <ev_viability>{ev_viability}</ev_viability>
    </environment>""")

        if cells:
            # write the cell agents
            for cell_idx in range(len(cells)):
                if isinstance(cells[cell_idx], elements.SecretoryCell):
                    of.write(cells[cell_idx].get_xml_description(target_source_points_per_cell))
                else:
                    of.write(cells[cell_idx].get_xml_description())

        # write the evs
        for ev in evs:
            ev.compute_msd(dt, degrees_of_freedom)
            # assign initial default velocity
            ev.vx = ev.velocity_ums * ev.base_direction_x
            ev.vy = ev.velocity_ums * ev.base_direction_y
            of.write(ev.get_xml_description(degrees_of_freedom, seconds_in_initial_state))

        # close the file
        of.write("""\n</states>""")

    print(f'file saved as: {f0}')

def parse_environment(tree):
    environment = {}
    for child in tree.xpath('/states/environment')[0]:
        tx = child.text
        if '.' in tx:
            environment[child.tag] = float(tx.replace('f', ''))
        else:
            environment[child.tag] = int(tx)
    return environment

def parse_agent(elem):
    data = {}
    try:
        for child in elem.getchildren():
            if child.text is None:
                value = None
            elif child.text == 'nan' or child.text == '-inf' or child.text == 'inf':
                value = None
                if 'errors_parsing' not in data:
                    data['errors_parsing'] = list()
                data['errors_parsing'].append(child.tag)
            elif '.' in child.text:
                if ',' in child.text:
                    value = [float(tx.replace('f', '')) for tx in child.text.split(',')]
                else:
                    value = float(child.text.replace('f', ''))
            elif child.tag == 'name':
                value = child.text
            else:
                value = int(child.text)
            data[child.tag] = value
        return data
    except ValueError as e:
        print('Error in',elem.xpath('child::name')[0].text,
            'id', elem.xpath('child::id')[0].text,'attribute',child.tag,
            'error:', e)
        return None
    except:
        print('Error in elem', elem.xpath('child::name')[0].text,
            'id', elem.xpath('child::id')[0].text,'child',child.tag,
            'error:', sys.exc_info()[0])
        return None

def parse_ev_agent(elem):
    data = {}
    try:
        for child in elem.getchildren():
            if child.text is None:
                value = None
            elif child.text == 'nan' or child.text == '-inf' or child.text == 'inf':
                value = None
                if 'errors_parsing' not in data:
                    data['errors_parsing'] = list()
                data['errors_parsing'].append(child.tag)
            elif '.' in child.text:
                if ',' in child.text:
                    value = [float(tx.replace('f', '')) for tx in child.text.split(',')]
                else:
                    value = float(child.text.replace('f', ''))
            elif child.tag == 'name':
                value = child.text
            else:
                value = int(child.text)
            if child.tag in ['id', 'x', 'y', 'z', 'age','name', 'radius_um']:
                data[child.tag] = value
        return data
    except ValueError as e:
        print('Error in',elem.xpath('child::name')[0].text,
            'id', elem.xpath('child::id')[0].text,'attribute',child.tag,
            'error:', e)
        return None
    except:
        print('Error in elem', elem.xpath('child::name')[0].text,
            'id', elem.xpath('child::id')[0].text,'child',child.tag,
            'error:', sys.exc_info()[0])
        return None

def parse_all_agents(tree):
    """
    evs[] contains a dictionary of tags:values as parsed from XML
    the grid based storages 'storage_name_g' contain the agents as Objects
    evs_with_errors contain the EVs as Objects
    """
    ciliary = {'objects':list(), 'dictionaries':list()}
    secretory = {'objects':list(), 'dictionaries':list()}
    evs = {'objects':list(), 'dictionaries':list()}
    evs_with_errors = {'objects':list(), 'dictionaries':list()}
    secretory_g = elements.GridStorage(5)
    ciliary_g = elements.GridStorage(5)
    all_cells_g = elements.GridStorage(5)
    evs_g = elements.GridStorage(5)

    for elem in tree.xpath('//xagent'):
        agent = parse_agent(elem)

        if agent['name'] == 'SecretoryCell':
            if 'source_points' in agent and agent['source_points'] > len(agent['source_points_xs']):
                agent['source_points_xs'] = agent['source_points_xs'][:agent['source_points']]
                agent['source_points_ys'] = agent['source_points_ys'][:agent['source_points']]
            del(agent['name'])
            secretory['dictionaries'].append(agent)
            extendedCell = elements.ExtendedCell(agent)
            secretory['objects'].append(extendedCell)
            secretory_g.store(extendedCell)
            all_cells_g.store(extendedCell)

        elif agent['name'] == 'CiliaryCell':
            del(agent['name'])
            ciliary['dictionaries'].append(agent)
            extendedCell = elements.ExtendedCell(agent)
            ciliary['objects'].append(extendedCell)
            ciliary_g.store(extendedCell)
            all_cells_g.store(extendedCell)

        elif agent['name'] == 'EV':
            ev = elements.EV()
            del(agent['name'])
            # for each agent variable parsed from XML, add the value read
            # to the EV object as an attribute. 
            # If NaNs were found, the agent is not added to the
            # general storage but to a list of evs with errors
            try:
                ev._id = int(agent['id'])
                for tag, value in agent.items():
                    setattr(ev, tag, value)

                if 'errors_parsing' in agent:
                    evs_with_errors['objects'].append(ev)
                    evs_with_errors['dictionaries'].append(agent)
                else:
                    evs_g.store(ev)
                    evs['objects'].append(ev)
                    evs['dictionaries'].append(agent)
            except AttributeError as err:
                print('Agent does not contain an attibute', agent)
                print(err)
    #print('Agents not parsed due to errors:',problems_count)
    return [secretory, ciliary, evs, evs_with_errors, secretory_g, ciliary_g,
            all_cells_g, evs_g]

def parse_ev_agents(filename):
    """
    evs contains a list of dictionary-based representations of the EVs as parsed from the XML
    """
    import xmltodict

    evs = []

    def parse_ev_agent_callback(_, element):
        if 'name' in element and element['name'] == 'EV':
            # add this EV to the collection
            # v3.6+ have individual defaultAt, disabledAt, disabledReason
            # older versions: becameAt[defaultAt, disabledAt, disabledReason]
            defaultAt = disabledAt = disabledReason = None
            if 'becameAt' in element:
                [defaultAt, disabledAt, 
                    disabledReason] = [s.strip() for s in element['becameAt'].split(',')]
            elif 'defaultAt' in element:
                defaultAt = element['defaultAt']
                disabledAt = element['disabledAt']
            #print(becameAt)
            evs.append({
                'id': int(element['id']),
                'x': float(element['x']),
                'y': float(element['y']),
                'radius_um': float(element['radius_um']),
                'age': float(element['age']),
                'defaultAt': int(defaultAt),
                'disabledAt': int(disabledAt)
            })
        return True

    print('Reading XML file in batches')
    with open(filename, 'rb') as xmlinput:
        xmltodict.parse(xmlinput, item_depth=2, item_callback=parse_ev_agent_callback)

    return evs

def parse_ev_agents_whole_file(filename):
    """
    evs contains a list of dictionary-based representations of the EVs as parsed from the XML
    """
    import xmltodict

    not_evs, evs = [], []
    print('Loading whole XML to memory')
    with open(filename, 'rb') as xmlinput:
        xd = xmltodict.parse(xmlinput.read())
    
    # remove 'environment' from xf['states']
    del(xd['states']['environment'])

    # remove not evs
    n = len(xd['states']['xagent'])
    for i in range(n):
        if xd['states']['xagent'][i]['name'] != 'EV':
            not_evs.append(i)
    for i in reversed(sorted(not_evs)):
        del xd['states']['xagent'][i]
    
    # 
    while(len(xd['states']['xagent']) > 1):
        o = xd['states']['xagent'].pop()
        defaultAt = disabledAt = disabledReason = None

        if 'becameAt' in o:
            [defaultAt, disabledAt, 
                disabledReason] = [s.strip() for s in o['becameAt'].split(',')]
        elif 'defaultAt' in o:
            defaultAt = o['defaultAt']
            disabledAt = o['disabledAt']
        evs.append({
                'id': int(o['id']),
                'x': float(o['x']),
                'y': float(o['y']),
                'radius_um': float(o['radius_um']),
                'age': float(o['age']),
                'defaultAt': int(defaultAt),
                'disabledAt': int(disabledAt)
            })
    evs.reverse()
    return evs

def read_xml(filename, ev_only=False, streaming=False):
    """
    if ev_only is enabled, returns a list of dictionaries and a list of objects
    both representing the EVs as parsed from XML plus another list which
    contains the evs producing errors while parsing
    """
    #filename = path.join(p, xml_file)

    if ev_only:
        if streaming:
            evs = parse_ev_agents(filename)
        else:
            evs = parse_ev_agents_whole_file(filename)
        return evs
    else:
        with open(filename,'r') as f:
            tree = etree.parse(f)
        itno = int(tree.xpath('/states/itno')[0].text)
        environment = parse_environment(tree)

        secretory, ciliary, evs, evs_with_errors, *grids = parse_all_agents(tree)
        print('read_xml len(grids):', len(grids))
        return [itno, environment, secretory, ciliary, evs, evs_with_errors, grids]

def agent_data_to_data_frame(agents, parameters=['id', 'x', 'y', 'age', 'radius_um']):
    columns = {}
    for parameter in parameters:
        columns[parameter] = list()

    for agent in agents:
        for parameter in parameters:
            if parameter in agent:
                columns[parameter].append(agent[parameter])
            else:
                columns[parameter].append(0)

    return pd.DataFrame.from_dict(columns)

def agent_data_to_data_frame_with_selection(agents, selected_elements, parameters=['id', 'x', 'y', 'age', 'radius_um']):
    df = agent_data_to_data_frame(agents, parameters)
    return [df, df.index[df.id.isin(selected_elements)].tolist()]

def agent_data_to_data_frame_without_selection(agents, selected_elements, parameters=['id', 'x', 'y', 'age', 'radius_um']):
    df = agent_data_to_data_frame(agents, parameters)
    df_oob = df[df.id.isin(selected_elements)]
    df_inbounds = df[~df.id.isin(selected_elements)]
    return [df_inbounds, df_oob]