import DelaunayPAN as DP
import matplotlib.pyplot as plt
import matplotlib.patches as pat

nodes = [ DP.Node(0.0, 0.0), DP.Node(-1.0, 2.0), DP.Node(0.0, 4.0), DP.Node(3.0, 5.0), DP.Node(5.0, 3.0), DP.Node(4.0, -1.0), DP.Node(2.0, 1.0) ]
boundaries = [ DP.Boundary([0, 1, 2, 3, 4, 5, 6], True) ]

nodes, elements = DP.delaunaymain(nodes, boundaries, 1.0, 100)

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)

for element in elements:
    ax.add_patch(pat.Polygon(xy = [(nodes[element.nodes[0]].x, nodes[element.nodes[0]].y), (nodes[element.nodes[1]].x, nodes[element.nodes[1]].y), (nodes[element.nodes[2]].x, nodes[element.nodes[2]].y)], fc = "blue", ec = "black"))

plt.xlim(-1.0, 5.0)
plt.ylim(-1.0, 5.0)
plt.axes().set_aspect('equal')
plt.show()