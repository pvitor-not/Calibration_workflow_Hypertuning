from collections import defaultdict
from lxml import etree as ET
from xml.dom import minidom
import h5py, io
import numpy as np


def etree_to_dict(tree):
    """
    Function to translate a xml-formatted content in the form of a dictionary

    text is accessible through the key '#text'
    attributes are accessible through the key '@name_of_the_attribute'
    """
    dico = {tree.tag: {} if tree.attrib else None}
    children = list(tree)

    if children:
        subdico = defaultdict(list)
        for child in map(etree_to_dict, children):
            for k, v in child.items():
                subdico[k].append(v)
        dico = {tree.tag: {k: v[0] if len(v) == 1 else v for k, v in subdico.items()}}
    if tree.attrib:
        dico[tree.tag].update(('@' + k, v) for k, v in tree.attrib.items())
    if tree.text:
        text = tree.text.strip()
        if children or tree.attrib:
            if text:
                dico[tree.tag]['#text'] = text
        else:
            dico[tree.tag] = text

    return dico


def dict_to_etree(dico: dict, parent=None, name=None):
    """
    Function to translate a dictionary in the form of a xml-formatted content

    text must be provided through the key '#text'
    attributes must be provided through the key '@name_of_the_attribute'

    : param dico:   xml content in dictionary format
    : param parent: parent element
    : return:       xml tree
    """

    if parent is None and name is None:
        if len(list(dico.keys())) != 1:
            raise Exception('There is more than one section in the xml: impossible to write this kind of file.')
        else:
            keys = list(dico.keys())
            parent = dict_to_etree(dico[keys[0]], name=keys[0])
            return parent

    parent = ET.Element(name)
    if isinstance(dico, dict):
        for key in dico.keys():
            if isinstance(dico[key], str):
                if key[0] == '@':
                    parent.set(key[1:], dico[key])
                elif key[0] == '#':
                    parent.text = dico[key]
                else:
                    child = ET.SubElement(parent, key)
                    child.text = dico[key]
            else:
                child = dict_to_etree(dico[key], parent=parent, name=key)
                if isinstance(child, list):
                    for c in child:
                        parent.append(c)
                else:
                    parent.append(child)
        return parent

    elif isinstance(dico, list):
        parent = []
        for ele in dico:
            if isinstance(ele, str):
                child = ET.Element(name)
                child.text = ele
                parent.append(child)
            else:
                parent += [dict_to_etree(ele, name=name)]

    elif isinstance(dico, str):
        child = ET.Element(name)
        child.text = dico
        parent.append(child)

    elif dico is None:
        parent = ET.Element(name)
        parent.text = ""

    else:
        raise Exception('Unknown situation found while building the xml file.')
    return parent


def xml_to_dict(file, encode='UTF-8'):
    """
    Function which returns a dictionary from an xml file

    Parameter:
        file: 	    path to the xml file
    Return:
       my_dict: 	dictionary with the xml structure
    """
    with io.open(file, 'r', encoding=encode) as xml:
        parser = ET.XMLParser(recover=True)
        xml_file = ET.fromstring(bytes(xml.read(), encoding=encode), parser)
        # xml_file = ET.fromstring(xml.read())
        my_dict = etree_to_dict(xml_file)

    return my_dict


def deep_keys(current_dict: dict, keys):
    """
    This is a recursive function to find the content of a deep dictionary

    Parameter:
        dict: initial dictionary
        keys: list of the successive keys to be followed
    Return:
        values: attribute of the last key
    """

    if len(keys) == 0:
        return current_dict
    else:
        return deep_keys(current_dict[keys[0]], keys[1:])


def reader_h5(file, data):
    with h5py.File(file, 'r') as f:
        data = [x for x in data.split('/') if x != '']
        results = deep_keys(f, data)
        return np.array(results)


def dict_to_xml_in_string(input_dict: dict, encode='UTF-8'):
    """
    Method to build the string equivalent of an xml file stored in a dictionary

    :param input_dict:  dictionary storing the xml file content
    :param encode:      encoding of the xml file

    return:             string corresponding to the content of the xml file
    """

    xml_in_string = minidom.parseString(ET.tostring(dict_to_etree(input_dict), encoding=encode))
    xml_in_string = xml_in_string.toprettyxml(encoding='UTF-8')
    xml_in_string = ''.join(map(chr, xml_in_string))

    return xml_in_string
