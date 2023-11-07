#!/usr/bin/python3

import json
from itertools import product

FLOW_FACTORY = {
    "vid": {"video_app": 1, "video_bitrate": [1280]},
    "if": {"internet_flow": 1, "if_bitrate": [12]},
    "bf": {"backlog_flow": 1},
}

DEFAULT_FLOW_FACTORY = {
    "vid": {"video_app": 0, "video_bitrate": []},
    "if": {"internet_flow": 0, "if_bitrate": []},
    "bf": {"backlog_flow": 0},
}

FLOW_TYPES = ["vid", "if", "bf"]

ALGORITHM_FACTORY = {
    "mt": (0, 0, 1, 0),
    "pf": (0, 0, 1, 1),
    "mlwdf": (1, 1, 1, 1),
}

ALGORITHM_PARAMS = ("algo_alpha", "algo_beta", "algo_epsilon", "algo_psi")

DEFAULT_UES_PER_SLICE = [5] * 10


def get_slice(flow_type, algorithm):
    if flow_type not in FLOW_TYPES:
        raise ValueError(f"Invalid flow type {flow_type}")
    if algorithm not in ALGORITHM_FACTORY:
        raise ValueError(f"Invalid algorithm {algorithm}")

    d_slice = {"n_slices": 10, "weight": 0.1}
    d_slice.update(FLOW_FACTORY[flow_type])
    for default_flow_type in FLOW_TYPES:
        if default_flow_type != flow_type:
            d_slice.update(DEFAULT_FLOW_FACTORY[default_flow_type])
    for k, v in zip(ALGORITHM_PARAMS, ALGORITHM_FACTORY[algorithm]):
        d_slice[k] = v

    return {"slices": [d_slice], "ues_per_slice": DEFAULT_UES_PER_SLICE}


if __name__ == "__main__":
    for flow_type, algorithm in product(FLOW_TYPES, ALGORITHM_FACTORY):
        with open(f"configs/{algorithm}-{flow_type}.json", "w", encoding="utf-8") as f:
            json.dump(get_slice(flow_type, algorithm), f, indent=2)
