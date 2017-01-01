# PCBTools
PCB Design Tools
This is a tool for Altium designer to automatically cut up power planes based on the location of power vias.
The idea is to create voronoi cells around each power via placed by the user and then merge adjacent voronoi cells
that connect to the same power net. Once all the adjacent voronoi cells have been merged, the user can assign the resultant
irregular polygon shapes to a particular power plane.
https://en.wikipedia.org/wiki/Voronoi_diagram
