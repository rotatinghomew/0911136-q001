"""村域矩形内的最近站分区面积权重（不依赖第三方库）。

站网疏密不均时，直接对站点取算术平均会让密集区重复计权。这里的做法是：
在村域矩形内按"离哪个站最近"划分面积，用各站所辖面积作权重，得到面积口径
的村域均值。只有一个站时，整村都归该站，权重为 1；此时并没有村内空间变化
信息，调用方必须把结果标成单站外推。
"""

import math

EPSILON = 1e-9


def rectangle_area(village):
    return (village["xmax"] - village["xmin"]) * (village["ymax"] - village["ymin"])


def point_in_rectangle(x, y, village):
    """半开区间 [xmin, xmax) × [ymin, ymax)，与村域不重叠的约定一致。"""
    return (village["xmin"] <= x < village["xmax"]
            and village["ymin"] <= y < village["ymax"])


def on_rectangle_boundary(x, y, village):
    return x in (village["xmin"], village["xmax"]) or y in (village["ymin"], village["ymax"])


def distance(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _intersection(p, q, a, b, c):
    denominator = a * (q[0] - p[0]) + b * (q[1] - p[1])
    if abs(denominator) < EPSILON:
        return p
    t = (c - a * p[0] - b * p[1]) / denominator
    return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))


def clip_halfplane(polygon, a, b, c):
    """Sutherland–Hodgman：保留凸多边形中 a*x + b*y <= c 的部分。"""
    if not polygon:
        return []
    output = []
    previous = polygon[-1]
    previous_inside = a * previous[0] + b * previous[1] <= c + EPSILON
    for current in polygon:
        current_inside = a * current[0] + b * current[1] <= c + EPSILON
        if current_inside:
            if not previous_inside:
                output.append(_intersection(previous, current, a, b, c))
            output.append(current)
        elif previous_inside:
            output.append(_intersection(previous, current, a, b, c))
        previous, previous_inside = current, current_inside
    return _deduplicate(output)


def _deduplicate(polygon):
    cleaned = []
    for point in polygon:
        if not cleaned or distance(cleaned[-1], point) > EPSILON:
            cleaned.append(point)
    if len(cleaned) > 1 and distance(cleaned[0], cleaned[-1]) <= EPSILON:
        cleaned.pop()
    return cleaned


def polygon_area(polygon):
    if len(polygon) < 3:
        return 0.0
    twice = 0.0
    for index, point in enumerate(polygon):
        following = polygon[(index + 1) % len(polygon)]
        twice += point[0] * following[1] - following[0] * point[1]
    return abs(twice) / 2.0


def thiessen_weights(points, village):
    """points 为村域内坐标；返回每个坐标在村域内所辖的面积（平方米）。

    坐标完全相同的站点先合并，否则它们的分区面积都会是 0；调用方负责先把
    重合站的雨量合并。
    """
    rectangle = [(village["xmin"], village["ymin"]), (village["xmax"], village["ymin"]),
                 (village["xmax"], village["ymax"]), (village["xmin"], village["ymax"])]
    if not points:
        return []
    if len(points) == 1:
        return [rectangle_area(village)]
    areas = []
    for index, own in enumerate(points):
        cell = list(rectangle)
        for other_index, other in enumerate(points):
            if index == other_index:
                continue
            normal = (other[0] - own[0], other[1] - own[1])
            midpoint = ((own[0] + other[0]) / 2.0, (own[1] + other[1]) / 2.0)
            cell = clip_halfplane(cell, normal[0], normal[1], normal[0] * midpoint[0] + normal[1] * midpoint[1])
            if not cell:
                break
        areas.append(polygon_area(cell))
    return areas


def max_corner_distance(points, village):
    """矩形四角到最近站点的最大距离；用来提示村内最远处离观测有多远。"""
    if not points:
        return None
    corners = [(village["xmin"], village["ymin"]), (village["xmax"], village["ymin"]),
               (village["xmax"], village["ymax"]), (village["xmin"], village["ymax"])]
    return max(min(distance(corner, point) for point in points) for corner in corners)


def nearest_neighbor_spacing(points):
    """站网典型最近邻间距的中位数；点太少时返回 None。"""
    coordinates = sorted({(round(point[0], 9), round(point[1], 9)) for point in points})
    if len(coordinates) < 2:
        return None
    nearest = []
    for index, point in enumerate(coordinates):
        others = [distance(point, other) for other_index, other in enumerate(coordinates)
                  if other_index != index]
        nearest.append(min(others))
    nearest.sort()
    middle = len(nearest) // 2
    if len(nearest) % 2:
        return nearest[middle]
    return (nearest[middle - 1] + nearest[middle]) / 2.0
