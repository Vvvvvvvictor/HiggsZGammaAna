import json
import os

import awkward
import numpy


_GLOBAL_TRACER = None


def set_global_tracer(tracer):
    global _GLOBAL_TRACER
    _GLOBAL_TRACER = tracer


def get_global_tracer():
    return _GLOBAL_TRACER


class YieldTracer:
    def __init__(self, path, truncate=False):
        self.path = os.path.abspath(path)
        self.weight_scale = 1.0

        if truncate:
            with open(self.path, "w", encoding="utf-8"):
                pass

    def set_weight_scale(self, weight_scale):
        self.weight_scale = float(weight_scale)

    @staticmethod
    def _to_numpy(values):
        if isinstance(values, awkward.Array):
            return awkward.to_numpy(values)
        return numpy.asarray(values)

    def record(self, stage, weights=None, mask=None, notes=""):
        if weights is None:
            if mask is None:
                raise ValueError("Either weights or mask must be provided for tracing.")
            if isinstance(mask, awkward.Array):
                n_events = len(mask)
            else:
                n_events = len(mask)
            weights_np = numpy.ones(n_events, dtype=numpy.float64)
        else:
            weights_np = self._to_numpy(weights).astype(numpy.float64, copy=False)

        if mask is None:
            mask_np = numpy.ones(len(weights_np), dtype=bool)
        else:
            if isinstance(mask, awkward.Array):
                mask_np = awkward.to_numpy(awkward.fill_none(mask, False))
            else:
                mask_np = numpy.asarray(mask)
            mask_np = mask_np.astype(bool, copy=False)

        if len(mask_np) != len(weights_np):
            raise ValueError(
                f"Mask and weight lengths differ: {len(mask_np)} vs {len(weights_np)} for stage {stage}."
            )

        selected = weights_np[mask_np]
        payload = {
            "stage": stage,
            "N_raw": int(numpy.sum(mask_np)),
            "sumw": float(numpy.sum(selected)),
            "sumw2": float(numpy.sum(numpy.square(selected))),
            "notes": notes,
        }

        with open(self.path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, ensure_ascii=True) + "\n")
