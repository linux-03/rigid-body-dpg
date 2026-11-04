"""Play a Recording in the browser.

The whole trajectory is sent to three.js as keyframes, so playback stays smooth
however long the simulation took, and the kernel is free while it plays.
"""
import numpy as np
import ipywidgets as widgets
import pythreejs as p3


def _quaternions(rotations):
    """(n, 3, 3) rotation matrices to continuous (n, 4) x, y, z, w quaternions."""
    out = np.empty((len(rotations), 4))
    for k, R in enumerate(rotations):
        K = np.array([
            [R[0, 0] - R[1, 1] - R[2, 2], 0, 0, 0],
            [R[0, 1] + R[1, 0], R[1, 1] - R[0, 0] - R[2, 2], 0, 0],
            [R[0, 2] + R[2, 0], R[1, 2] + R[2, 1], R[2, 2] - R[0, 0] - R[1, 1], 0],
            [R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1], R[0, 0] + R[1, 1] + R[2, 2]],
        ]) / 3.0
        q = np.linalg.eigh(K)[1][:, -1]
        if k and q @ out[k - 1] < 0:
            q = -q
        out[k] = q
    return out


def _floats(values):
    return np.asarray(values, dtype=np.float32).ravel()


def _rod_tracks(name, ends, times):
    start = np.array([e[0] for e in ends])
    stop = np.array([e[1] for e in ends])
    axis = stop - start
    length = np.linalg.norm(axis, axis=1)
    y = axis / length[:, None]
    # Shortest rotation taking the cylinder's +y axis onto the rod.
    quats = np.column_stack([y[:, 2], np.zeros(len(y)), -y[:, 0], 1.0 + y[:, 1]])
    quats /= np.linalg.norm(quats, axis=1, keepdims=True)
    scale = np.column_stack([np.ones(len(y)), length, np.ones(len(y))])
    return [
        p3.VectorKeyframeTrack(name=f"{name}.position", times=times, values=_floats(0.5 * (start + stop))),
        p3.QuaternionKeyframeTrack(name=f"{name}.quaternion", times=times, values=_floats(quats)),
        p3.VectorKeyframeTrack(name=f"{name}.scale", times=times, values=_floats(scale)),
    ]


def _extent(rec):
    """Centre and diameter of a sphere containing the whole motion."""
    points = [rec.positions.reshape(-1, 3)]
    for body_index, body in enumerate(rec.world.bodies):
        radius = np.max(np.linalg.norm(body.vertices, axis=1))
        centers = rec.positions[:, body_index]
        for axis in np.eye(3):
            points += [centers + radius * axis, centers - radius * axis]
    for frames in (rec.link_ends, rec.spring_ends):
        for frame in frames:
            if frame:
                points.append(np.array(frame).reshape(-1, 3))
    points = np.concatenate(points)
    center = 0.5 * (points.min(axis=0) + points.max(axis=0))
    return center, max(2.0 * float(np.max(np.linalg.norm(points - center, axis=1))), 0.5)


def show(rec, speed=1.0, width=760, height=460, background="#f6f5f1", autoplay=True):
    """Build the animated scene for a Recording."""
    world = rec.world
    times = np.asarray(rec.times / speed, dtype=np.float32)
    objects, tracks = [], []

    for i, body in enumerate(world.bodies):
        geometry = p3.BufferGeometry(attributes={
            "position": p3.BufferAttribute(array=np.asarray(body.vertices, dtype=np.float32)),
            "normal": p3.BufferAttribute(array=np.asarray(body.normals, dtype=np.float32)),
        })
        material = p3.MeshStandardMaterial(color=body.color, roughness=0.55, metalness=0.05, side="DoubleSide")
        objects.append(p3.Mesh(geometry, material, name=body.name))
        tracks += [
            p3.VectorKeyframeTrack(name=f"{body.name}.position", times=times, values=_floats(rec.positions[:, i])),
            p3.QuaternionKeyframeTrack(name=f"{body.name}.quaternion", times=times,
                                       values=_floats(_quaternions(rec.rotations[:, i]))),
        ]

    center, size = _extent(rec)
    rod_radius = 0.008 * size
    rod_geometry = p3.CylinderBufferGeometry(radiusTop=rod_radius, radiusBottom=rod_radius, height=1, radialSegments=12)
    rod_material = p3.MeshStandardMaterial(color="#3a3a3a", roughness=0.4)
    for k in range(len(world.links)):
        name = f"rod{k}"
        objects.append(p3.Mesh(rod_geometry, rod_material, name=name))
        tracks += _rod_tracks(name, [frame[k] for frame in rec.link_ends], times)

    # A helix of unit length along y; scaling it in y stretches the coils like a real spring.
    turns = 10
    angle = np.linspace(0, 2 * np.pi * turns, 40 * turns)
    coil_radius = 0.025 * size
    helix = np.column_stack([coil_radius * np.cos(angle), np.linspace(-0.5, 0.5, len(angle)), coil_radius * np.sin(angle)])
    spring_geometry = p3.BufferGeometry(attributes={"position": p3.BufferAttribute(array=np.asarray(helix, dtype=np.float32))})
    spring_material = p3.LineBasicMaterial(color="#30638e")
    for k in range(len(world.springs)):
        name = f"spring{k}"
        objects.append(p3.Line(spring_geometry, spring_material, name=name))
        tracks += _rod_tracks(name, [frame[k] for frame in rec.spring_ends], times)

    group = p3.Group(children=objects)
    clip = p3.AnimationClip(tracks=tracks, duration=float(times[-1]))
    action = p3.AnimationAction(p3.AnimationMixer(group), clip, group)

    up = -world.gravity / np.linalg.norm(world.gravity) if np.any(world.gravity) else np.array([0.0, 0.0, 1.0])
    side = np.cross(up, [1.0, 0.0, 0.0])
    side = side if np.linalg.norm(side) > 1e-6 else np.cross(up, [0.0, 1.0, 0.0])
    side /= np.linalg.norm(side)
    view = side + 0.35 * up + 0.6 * np.cross(side, up)
    fov = 35.0
    # Far enough that the sphere around the motion fits the view's smaller (vertical) extent.
    eye = center + (0.5 * size / np.sin(np.radians(0.5 * fov))) * view / np.linalg.norm(view)
    camera = p3.PerspectiveCamera(
        position=tuple(eye), up=tuple(up), aspect=width / height, fov=fov, near=0.01 * size, far=50 * size,
        children=[p3.DirectionalLight(color="white", intensity=0.75, position=(3, 5, 4))],
    )
    grid = p3.GridHelper(4 * size, 20, "#b9b6ab", "#dcd9cf")
    if abs(up[1]) < 0.5:  # GridHelper lies in the xz plane
        grid.quaternion = (float(np.sqrt(0.5)), 0.0, 0.0, float(np.sqrt(0.5)))
    grid.position = tuple(center - 0.5 * size * up)
    scene = p3.Scene(children=[camera, group, grid, p3.AmbientLight(color="white", intensity=0.55)], background=background)
    renderer = p3.Renderer(
        camera=camera, scene=scene, width=width, height=height, antialias=True,
        controls=[p3.OrbitControls(controlling=camera, target=tuple(center))],
    )
    hint = widgets.HTML("<small>Drag to rotate, scroll to zoom. Use the buttons to play, pause and stop.</small>")
    out = widgets.VBox([renderer, widgets.HBox([action, hint])])
    if autoplay:
        action.play()
    return out
