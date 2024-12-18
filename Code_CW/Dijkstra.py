import numpy as np


class Node_Distance:

    def __init__(self, name, dist):
        self.name = name
        self.dist = dist


class Graph:

    def __init__(self, list_OFS, node_count, nb_data, nb_cells):
        self.nb_cells = nb_cells
        self.nb_data = nb_data
        self.node_count = node_count
        self.OFS = list_OFS

    def Dijkstras_Shortest_Path(self, source, target):

        # Initialize the distance of all the nodes from source to infinity
        distance = [999999999999] * self.node_count
        # Distance of source node to itself is 0
        distance[source] = 0

        # Create a dictionary of { node, distance_from_source }
        dict_node_length = {}

        # Initialization with the first cell
        for j, length_to_adjnode in enumerate(self.OFS[0][0]):
            distance[j] = length_to_adjnode
            dict_node_length[j] = length_to_adjnode

        while dict_node_length:

            # Get the key for the smallest value in the dictionary
            # i.e Get the node with the shortest distance from the source
            source_node = min(dict_node_length, key=lambda k: dict_node_length[k])
            if source_node == target:
                return distance[target]
            del dict_node_length[source_node]

            cell_id = source_node // (self.nb_data + 1)
            top_marker = source_node % (self.nb_data + 1)

            # The destination of the last cell is necessarily the last node
            if cell_id == self.nb_cells - 2:
                length_to_adjnode = self.OFS[-1][top_marker][-1]
                adjnode = (self.nb_cells - 1) * (self.nb_data + 1)

                # Edge relaxation
                if distance[adjnode] > distance[source_node] + length_to_adjnode:
                    distance[adjnode] = distance[source_node] + length_to_adjnode
                    dict_node_length[adjnode] = distance[adjnode]

            # for a cell in the middle in the sequence
            else:
                for j, length_to_adjnode in enumerate(self.OFS[cell_id + 1][top_marker]):
                    adjnode = (cell_id + 1) * (self.nb_data + 1) + top_marker + j

                    # Edge relaxation
                    if distance[adjnode] > distance[source_node] + length_to_adjnode:
                        distance[adjnode] = distance[source_node] + length_to_adjnode
                        dict_node_length[adjnode] = distance[adjnode]

        return distance[target]
