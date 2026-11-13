"""处理经过质量标记的站点观测。"""

from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True)
class Reading:
    station_id: str
    village_id: str
    observed_at: datetime
    millimeters: float | None
    quality: str = "valid"


def valid_readings(readings: list[Reading]) -> list[Reading]:
    return [item for item in readings if item.quality == "valid" and item.millimeters is not None]
