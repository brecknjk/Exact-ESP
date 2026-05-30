import matplotlib.pyplot as plt

def read_gl_file(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Read counts
    num_nodes = int(lines[0])
    num_edges = int(lines[1])
    
    # Read node coordinates
    nodes = []
    for i in range(2, 2 + num_nodes):
        x, y = map(float, lines[i].split())
        nodes.append((x, y))
    
    # Read edges
    edges = []
    for i in range(2 + num_nodes, 2 + num_nodes + num_edges):
        n1, n2, n3, n4 = map(int, lines[i].split())
        edges.append((n1, n2))
    
    return nodes, edges

def visualize_graph(nodes, edges):
    plt.figure(figsize=(6, 6))
    
    # Plot edges
    for (n1, n2) in edges:
        x_values = [nodes[n1][0], nodes[n2][0]]
        y_values = [nodes[n1][1], nodes[n2][1]]
        plt.plot(x_values, y_values, 'b-', linewidth=1)
    
    # Plot nodes
    xs, ys = zip(*nodes)
    plt.scatter(xs, ys, c='red', s=50, zorder=5)
    
    # Label nodes
    for i, (x, y) in enumerate(nodes):
        plt.text(x, y, f'{i}', fontsize=9, ha='right', va='bottom')
    
    plt.axis('equal')
    plt.title('Graph Visualization')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

def main():
    filename = input("Enter the path to the .gl file: ").strip()
    nodes, edges = read_gl_file(filename)
    visualize_graph(nodes, edges)

if __name__ == '__main__':
    main()
