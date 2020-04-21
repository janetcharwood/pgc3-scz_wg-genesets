#!/usr/bin/env python
def make_element_dict(elem_type, attr_dict, ontology, elem_dict):
    for k, v in attr_dict.items():
        if ontology in v.get('namespace'):
            # construct element
            el = Element(elem_type)
            # print v.get('namespace')
            el_id = ''.join(v.get('id'))
            # print el_id
            # set this as id in Element object and attribute dictionary for term
            el.set('id', el_id)
            # print out  the element type and its attibutes
            # print '', el.tag,el.attrib
            # store Element object for term in elem_dict
            # bp_elem_dict[el_id] = el
            elem_dict[el_id] = el

##########################################################################
