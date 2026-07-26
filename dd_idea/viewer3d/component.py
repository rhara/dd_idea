"""Double-buffered Streamlit embedding for py3Dmol scenes, vendored from
`dd_viewer.component` (renamed component id only, to keep it distinct from
the sibling `dd_seqalign` app's own `plviewer_3d` instance if both ever run
in the same browser session).

`st.iframe(html)` replaces the iframe's srcdoc on every Streamlit rerun; the
browser tears down the old document immediately and the new one takes a
moment to load 3Dmol.js, re-parse the model, and re-render, producing a
visible blank flash on every widget interaction (a residue-table selection,
a style checkbox, ...).

`view3d` embeds scenes through a small static Streamlit component
(`frontend/index.html`) instead. Because it's a *static* component, its own
outer iframe has a fixed local URL and is loaded exactly once; Streamlit
talks to an already-running instance via postMessage rather than reloading
it, so the component's JS can stay alive across reruns. It keeps two hidden/
visible inner iframes and only swaps which one is visible once the new scene
reports (via the postMessage `htmlpatch.html_with_camera_events` injects)
that it has actually finished painting -- so updates read as a quick cross-
fade instead of a flash to blank. That same always-alive JS context is also
what lets the camera position (rotation/zoom) survive across scenes: it's
kept in a plain JS variable in the component, not in the scene's own (short-
lived, storage-partitioned) iframe.
"""
import os

import streamlit.components.v1 as components

_FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend")
_component = components.declare_component("dd_idea_3d", path=_FRONTEND_DIR)


def view3d(html: str, height: int = 650, key: str = "dd_idea_3d_view", reset_camera_token: int = 0):
    """Embed a py3Dmol scene's HTML (typically already passed through
    `htmlpatch.html_with_camera_events`) via the double-buffered component.

    `key` should stay constant across reruns (the default already does) --
    it's what lets Streamlit recognize this as the same live component
    instance rather than tearing down and recreating it, which would defeat
    the whole point. `reset_camera_token` forgets the saved camera position
    (falling back to the scene's own zoomTo fit) whenever it changes from
    the previous call -- wire it to e.g. a counter incremented by a "reset
    view" button.

    Returns the latest reported camera view (a `getView()`-shaped list, or
    `None` before the component has ever reported one) -- pass it back into
    the *next* call's `html` via `htmlpatch.html_with_initial_view` before
    calling this again. This isn't optional belt-and-suspenders: a live
    comparison against a real multi-protein report showed Streamlit does
    sometimes tear down and recreate this component's own outer iframe
    (and with it, the in-memory `lastView`/`activeKey` its `frontend/
    index.html` keeps alive across reruns to restore the camera) on an
    ordinary widget-triggered rerun -- confirmed by a marker stashed on the
    iframe's `contentWindow` disappearing across the rerun, and the next
    scene coming back at its own default `zoomTo()` fit instead of the
    camera position just set before that rerun. Feeding the view back in
    from `st.session_state` through the HTML itself survives that
    regardless of what the component's own JS memory does.
    """
    return _component(html=html, height=height, key=key, reset_token=reset_camera_token, default=None)
