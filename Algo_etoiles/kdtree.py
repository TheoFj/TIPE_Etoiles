import heapq

class KDtree:
    def __init__(self, star, left, right):
        self.star = star
        self.left = left
        self.right = right

def build_tree(lst, dim, dir): #O(nlog(n)^2)
    if lst==[]:
        return None
    lst.sort(key = lambda star: star.pos[dir])
    mid = len(lst)//2
    medianstar = lst[mid]
    lows, highs = lst[:mid], lst[(mid+1):]
    return KDtree(medianstar, build_tree(lows, dim, (dir+1)%dim), build_tree(highs, dim, (dir+1)%dim))

def distsqr(p1, p2, dim):
    return sum([(p2[i]-p1[i])**2 for i in range(dim)])

def nearest_nstars(node, targetpos, n, dim, dir, heap=[], withtarget=False):
    if node is None:
        return
    
    dist = distsqr(node.star.pos, targetpos, dim)
    #On met des distances negaties pour simuler un max-heap : on veut que le pop retire à chaque fois celui dont la distance à targetpos est maximale càd la moins bonne etoile
    if (dist!=0 or withtarget):
        heapq.heappush(heap, (-dist, node.star.pos, node.star)) #le star.pos sert a differencier des etoiles pouvant etre a la meme distance de target car la heapq comparant en ordre lexicographique, cela evite de comparer les objets etoile
    if len(heap) > n:
        heapq.heappop(heap)

    diff = targetpos[dir] - node.star.pos[dir]
    primary_branch = node.left if diff < 0 else node.right
    secondary_branch = node.right if diff < 0 else node.left
    
    # Explorer d'abord la branche principale
    nearest_nstars(primary_branch, targetpos, n, dim, (dir + 1)%dim, heap)

    # Explorer l'autre branche si nécessaire
    if len(heap) < n or (diff)**2 <= -heap[0][0]:
        nearest_nstars(secondary_branch, targetpos, n, dim, (dir + 1)%dim, heap)

    return heap