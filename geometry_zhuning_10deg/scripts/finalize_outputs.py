# -*- coding: utf-8 -*-
"""整理真实 Abaqus 导图：只裁去空白边框、加文字，不改变几何轮廓。"""
from pathlib import Path
import json
import math
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent.parent  # 脚本位于 scripts\，数据根目录为上级
latest = json.loads((ROOT / 'latest_run.json').read_text(encoding='utf-8'))
run = Path(latest['run_directory'])
render = json.loads((run / 'render_status.json').read_text(encoding='utf-8'))
assert render['status'] == 'COMPLETE'
font_path = r'C:\Windows\Fonts\msyh.ttc'
font = ImageFont.truetype(font_path, 32)
small = ImageFont.truetype(font_path, 24)
results = []
for source, name, title in [
    ('abaqus_section_full.png', 'model_section_full.png', '10°扇区子午截面 · Abaqus/CAE 实体渲染'),
    ('abaqus_section_joint.png', 'model_section_joint.png', '阀座—内衬连接区 · Abaqus/CAE 实体渲染'),
    ('abaqus_sector_3d.png', 'model_sector_3d.png', '10°基体实体 · Abaqus/CAE 实体渲染')]:
    original = Image.open(run / source).convert('RGB')
    a = np.asarray(original).astype(float)
    r, g, b = a[:, :, 0], a[:, :, 1], a[:, :, 2]
    blue = (r < 180) & (g > 70) & (g < 210) & (b > 95) & (b > r * 1.2)
    gold = (r > 120) & (g > 70) & (b < 145) & (r > g * 1.08) & (g > b * 1.15)
    ys, xs = np.where(blue | gold)
    assert len(xs) > 1000
    bounds = (max(0, int(xs.min()) - 20), max(0, int(ys.min()) - 20),
              min(original.width, int(xs.max()) + 21), min(original.height, int(ys.max()) + 21))
    crop = original.crop(bounds)
    canvas = Image.new('RGB', (max(crop.width + 60, 1160), crop.height + 180), 'white')
    draw = ImageDraw.Draw(canvas)
    draw.text((30, 18), title, fill='#202B37', font=font)
    canvas.paste(crop, ((canvas.width - crop.width) // 2, 78))
    y = canvas.height - 66
    draw.rectangle((30, y + 4, 55, y + 26), fill='#398FBD')
    draw.text((68, y), 'PA6 内衬', fill='#202B37', font=small)
    draw.rectangle((255, y + 4, 280, y + 26), fill='#D9A14C')
    draw.text((294, y), '6061-T6 铝阀座', fill='#202B37', font=small)
    draw.text((560, y), '未生成复材层；未标注接头尺寸为拟合值；非 AI 图像', fill='#475569', font=small)
    canvas.save(ROOT / name)
    results.append({'image': str(ROOT / name), 'raw_abaqus_image': str(run / source),
                    'operations': ['裁去空白及坐标图标', '添加标题与材料图例'],
                    'geometry_pixels_modified': False, 'ai_generated': False})

# 对照78 L公称值，计算两端瓶口端平面之间的封闭包络减去实体材料体积。
profiles = json.loads((ROOT / 'profiles.json').read_text(encoding='utf-8'))
validation = json.loads((run / 'geometry_validation.json').read_text(encoding='utf-8'))
d = profiles['parameters']['paper_dimensions']
a, b, h = profiles['derived']['ellipse_a'], profiles['derived']['ellipse_b'], profiles['derived']['tangent_y']
fc = profiles['derived']['fillet_center']
py = profiles['derived']['ellipse_fillet_point'][1]
rf = profiles['parameters']['fitted_joint_dimensions']['neck_outer_fillet_radius']
neck, end = d['pole_outer_diameter'] / 2, d['overall_length'] / 2
outer_half_integral = a * a * h
outer_half_integral += quad(lambda y: a * a * (1 - ((y - h) / b) ** 2), h, py)[0]
outer_half_integral += quad(lambda y: (fc[0] - math.sqrt(max(0, rf * rf - (y - fc[1]) ** 2))) ** 2, py, fc[1])[0]
outer_half_integral += neck * neck * (end - fc[1])
outer_volume = 2 * math.pi * outer_half_integral
material_volume = validation['merged_volume_mm3'] * 360 / d['sector_angle_deg']
capacity = (outer_volume - material_volume) / 1e6
capacity_report = {'geometric_cavity_litre': capacity, 'paper_nominal_litre': 78.0,
    'difference_litre': capacity - 78.0, 'difference_percent': (capacity / 78 - 1) * 100,
    'definition': '两端瓶口平面封闭后的基体内部容积；未扣除实际瓶阀/盲堵占据体积。',
    'status': 'NOT_A_NOMINAL_CAPACITY_CERTIFICATION',
    'note': '论文仅给部分接头尺寸，本模型优先满足已标注尺寸。拟合接头侵入容腔使实际几何容积与78 L公称值存在差异；因此不能宣称已完全复现论文或满足国标公称容积要求。'}
(run / 'capacity_check.json').write_text(json.dumps(capacity_report, ensure_ascii=False, indent=2), encoding='utf-8')
manifest = {'cae_file': latest['cae_file'], 'sat_file': str(run / 'zhuning_base_10deg.sat'),
    'validation_file': str(run / 'geometry_validation.json'), 'capacity_check': capacity_report,
    'images': results, 'geometry_status': validation['status'],
    'reproduction_scope': '已标注尺寸一致；接头缺失尺寸为有记录的近似，非完全同一模型。'}
(ROOT / 'delivery_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'images': [q['image'] for q in results], 'capacity_check': capacity_report}, ensure_ascii=False, indent=2))
