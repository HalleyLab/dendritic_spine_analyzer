# 用 Fiji/ImageJ 手工复现树突棘阳性面积／树突长度

全程使用菜单和手工 ROI，不需要宏。目标是二维荧光阳性棘候选面积 A，除以树突中心线弧长 L，单位 µm²/µm。它不是棘数量、三维体积或细胞膜表面积。

建议先按第 1–14 步独立处理 CZI，再按文末的“核对程序结果”复算。前者检查分割是否合理，后者检查比例尺、像素计数与除法是否一致。

## 1. 打开原始 CZI，核对要处理的图像

1. 在 Fiji 选择 `Plugins > Bio-Formats > Bio-Formats Importer`，打开 `your_image.czi`。
2. 通道、场景、时间点必须与 Notebook 输入一致。当前文件读取到一个 AF546 通道、1024×1024 的二维图像。
3. 如果换成带 Z 层的文件，在目标通道/时间点上使用 `Image > Stacks > Z Project...`，选择与程序相同的层范围和 `Max Intensity`。程序的 Z 索引从 0 开始，ImageJ 的层号通常从 1 开始。
4. 当前 final.czi 已是二维数据，若没有 Z 堆栈就直接使用，不要再次人为“投影”。用 `Image > Duplicate...` 建立工作副本，保留原始 CZI。

核验：原图的树突走向、胞体位置应与程序左侧原图相同。[Bio-Formats 导入说明](https://bio-formats.readthedocs.io/en/stable/users/imagej/load-images.html)

## 2. 核对比例尺和原始位深

1. 打开 `Image > Properties...`，检查 `Pixel width`、`Pixel height` 和 `Unit of length`。
2. 当前文件 X、Y 均为 **0.0990358382936508 µm/pixel**，单位可显示为 `um`、`micron` 或 `µm`。实际以你当前 CZI 元数据为准。
3. 保持原始 16-bit 图像用于阈值分析。显示为红色只是 LUT，不代表图像已经变成 RGB。
4. 不要拿截图大小推算比例尺；不要开启影响所有图像的 Global scale。若标定不符，先解决导入问题，再测量。

核验：本文件每个原始像素面积为约 **0.0098080973 µm²**。程序面积采用 X×Y，即使 X/Y 不相等也分别使用。[图像属性说明](https://imagej.net/ij/docs/guide/146-28.html)

## 3. 选靠近胞体分叉的一条树突，并记录真实长度

1. 在原图找到胞体及其附近的主树突分叉，沿其中一条分支向远离胞体的方向选段。避开胞体本身、分叉膨大区域和交叉纤维；不要跨到另一条分支。
2. 选择工具栏的 `Segmented Line`。由近端向远端沿主干中心依次点选，弯曲处加点，双击结束。
3. 起点偏移必须与程序定义一致：程序默认从**你点选的首个中心线点**向远端偏移 2 µm，然后截取 15 µm；它没有自动识别胞体或分叉。手工操作时，先确定偏移后的正式测量起点，再画正式中心线。
4. 在 `Analyze > Set Measurements...` 关闭 `Limit to threshold`，设置足够小数位；使用 `Analyze > Measure` 读取 `Length`。拖动末端控制点，使 L 接近 15 µm。
5. 手工画线存在定位误差，可事先约定例如一个像素尺寸的长度容差。最终分母使用**实际测得的 L**，不要把 14.8 µm 强行当成 15 µm。
6. 不做 `Fit Spline`，以便与程序的折线弧长定义对照。按 T 加入 `Analyze > Tools > ROI Manager...`，重命名为 `Centerline`，保存 ROI。

核验：线沿轴向穿过主干中央；测量的是弯曲路径长度，不是两个端点的直线距离。[选择与折线操作](https://imagej.net/ij/docs/guide/146-27.html)

## 4. 横向展直只用于观察

1. 另复制一份工作图，恢复刚才的中心线 ROI。
2. 双击线工具设置线宽，覆盖主干和两侧棘。若参考程序默认每侧 6 µm，本文件总宽约为 `12 / 0.099 ≈ 121 pixels`。
3. 使用 `Edit > Selection > Straighten...` 得到水平图，另存观察图。
4. 返回原图进行后续正式面积测量。ImageJ 展直可能拟合样条并重采样，不能要求展直图的像素面积与原图严格相同。

核验：横图只是帮助检查主干边界，最终 A 始终来自原始坐标的阳性掩膜。[Straighten 说明](https://imagej.net/ij/docs/guide/146-27.html)

## 5. 建立固定长度的分析范围

1. 在原图用 Polygon/Freehand ROI 勾画中心线两侧的带状范围，两端应对应正式测量的起止位置，并大致垂直于当地树突方向。
2. 将该 ROI 加入管理器，命名 `Corridor`。选区应覆盖完整棘头、棘颈，但尽量避开邻近树突。
3. 在 `Process > Binary > Options...` 勾选 `Black background`，随后 `Edit > Selection > Create Mask`，得到范围掩膜，白色区域数值为 255、黑色为 0。
4. 保存 `Corridor.tif`。以后所有 AND 运算都使用相同尺寸、相同原点的图像。

核验：矩形框只是定位辅助，实际范围还要限制轴向起止位置。**不要直接把程序裁剪尺寸的 `03_corridor.tif` 与 1024×1024 原图做 AND。**严格复算时，`03_corridor.tif` 只能与同尺寸的 `01_raw_crop.tif`／`02_corrected_crop.tif` 配套；独立流程则自行制作全尺寸范围 Mask。[二值掩膜说明](https://imagej.net/ij/docs/guide/146-29.html)

## 6. 调荧光阈值，保留棘头和棘颈

1. 第一次对照建议将程序的 Gaussian sigma 和 Background radius 均设为 0，ImageJ 也先不加滤波或减背景。
2. 在原始位深的工作副本上打开 `Image > Adjust > Threshold...`，点击 `Set`，下限填程序的 **Spine fluorescence threshold**，16-bit 上限可设 65535。
3. 放大观察：真实棘的较暗外缘和细颈若没有被阈值高亮，逐步降低下限。例如可把 32000、24000、18000、12000 作为本图的比较点，**这些只是试验值，不是推荐的统一正确阈值**。
4. 同时看背景和旁边纤维。降低阈值会纳入更多暗信号，也可能把背景连接成假结构。不要以“面积越大越好”选择阈值。
5. 用 `Process > Binary > Convert to Mask` 转为二值图。确认白色=255、黑色=0；必要时查看鼠标所在像素的数值，不能只看 LUT 的颜色。
6. 用 `Process > Image Calculator...` 将阳性 Mask 与 Corridor 做 `AND`，勾选 `Create new window`，得到范围内阳性 Mask，保存为 `Positive.tif`。

核验：阈值使用原始灰度，不是显示亮度。不要点击亮度/对比度的 Apply 来修改数据。若程序启用了预处理，严格比阈值时使用同尺寸的 `02_corrected_crop.tif`；程序的矩形开运算背景与 Fiji rolling ball 不等价。[Threshold 与显示设置](https://imagej.net/ij/docs/guide/146-28.html)

## 7. 独立勾画主干，防止把棘吞进主干

1. 在原图放大该段，用 Polygon 或 Freehand Selection 沿**连续主干的真实两侧边缘**勾画，不要绕到棘头外缘。
2. 加入 ROI Manager，命名 `ShaftDomain`，保存 ROI；使用 `Create Mask` 得到主干几何 Mask。
3. 用阳性 Mask AND ShaftDomain，得到主干阳性 Mask，命名 `ShaftPositive`。
4. 主干与棘都发荧光，ImageJ 单一阈值不会自动判断解剖归属；这一步需要结合形态人工确定。
5. 新版程序将 `Shaft reference` 与 `Spine fluorescence threshold` 分开。比较不同棘阈值时，在 ImageJ 保持同一主干 ROI，也能避免阈值变化导致主干同时变宽。

核验：若原本认为是棘的部分落在绿色主干范围内，应先改主干边界。仅降低棘阈值不能把已经划入主干的区域变回棘。

## 8. 去掉干扰，并检查连接关系

1. 保留 Positive 原副本，在另一个副本上处理。
2. 对不属于这条树突的纤维、侧枝、明显噪点，用 ROI 选中，确认背景色为黑色后使用 `Edit > Clear`，使对应像素变为 0。
3. 仅保留与该主干相连的阳性结构。孤立亮点先检查是否有较暗棘颈，不直接全部计为棘。
4. 有原始 Z 堆栈时，回到单层核实投影交叉；只有二维最终图时，无法可靠分开所有轴向重叠。
5. 不要误删主干或棘颈，使真实棘断开；保存剔除区域 ROI，以便复核。

核验：程序采用 8 邻域连接。你手工认定的归属规则应记录下来，尤其是暗颈与交叉结构。

## 9. 扣除主干，得到筛选前的棘候选

1. 复制 ShaftDomain Mask，使用 `Edit > Invert` 数值反相：原主干内部变成 0，外部变成 255。
2. 注意不是 `Invert LUT`；后者只换显示颜色，不改变 0/255 数值。
3. 用清理后的阳性 Mask AND 反相主干 Mask，得到 `CandidatesBeforeFilters.tif`。
4. 不做 Fill Holes，不把棘轮廓包住的所有空白面积填进去。目标是荧光阳性像素面积。

核验：这张图应已去掉主干，但可能仍含过小噪声、过大的连通块或较长侧枝。[二值与图像运算](https://imagej.net/ij/docs/guide/146-29.html)

## 10. 应用面积和延伸限制，并记录排除原因

1. 在候选 Mask 上使用 `Analyze > Analyze Particles...`，Size 输入与程序一致的范围。当前默认最小 0.001 µm²、最大面积 0 表示不设上限，因此 ImageJ 的 Size 应输入 `0.001-Infinity`，不要把最大面积填成 0。
2. 图像标定为微米后，应使用平方微米，避免勾选 Pixel Units。Circularity 设 0–1，不额外按圆度筛掉细长棘。
3. 不勾 Include Holes；程序不会仅因碰到选区边缘就自动全部排除，所以对照时不要擅自勾 Exclude on Edges。使用 Outlines/ROI Manager 辅助定位对象。
4. 保留原候选 Mask 的像素，不直接拿填充后的轮廓图替代它。将确认应排除的对象在候选 Mask 副本上整块清为 0，保留孔洞。
5. 延伸上限默认 10 µm，是从主干边缘向外的最远距离。用局部法向线辅助测量；超过上限时，程序排除整个突出部连通区域，不是截掉超出部分后保留余下侧枝。
6. 相邻多个棘若连成一片，可能被合并为一个大区域。结合原图检查，不能把这种连通区域数当成棘数；放宽面积/延伸限值后也要排除真正侧枝。

核验：若 CandidatesBeforeFilters 很多而最终很少，原因在筛选或连接关系，不能只继续降阈值。新版程序的 `Before area/reach filters`、`Rejected / detached signal` 视图和面积统计用于定位这一问题。[粒子分析说明](https://imagej.net/ij/docs/guide/146-30.html)

## 11. 测最终棘阳性面积 A

1. 在最终棘 Mask 上使用 `Edit > Selection > Select None`，避免只测某个残留 ROI。
2. 在 Threshold 窗口设置 **255–255**。
3. 打开 `Analyze > Set Measurements...`，勾选 **Area、Limit to threshold**，Redirect to 设 **None**。
4. `Analyze > Measure`，记录结果表中 `Area` 为 A，单位 µm²，保存结果表。

核验：不勾 Limit to threshold 可能测到整张矩形面积；%Area 是比例，不是 A；彩色标记 PNG 也不能代替二值 Mask 测面积。[测量说明](https://imagej.net/ij/docs/guide/146-30.html)

## 12. 再测长度 L，并计算 A/L

1. 回到原始灰度图，恢复 Centerline ROI。
2. 关闭 Limit to threshold；Measure 读取 `Length`，确认原图标定没有变化。
3. 用实际测得的 A 除以实际 L，得到 µm²/µm。保存 A、L、比值、阈值和筛选规则，而不是只保存最终比值。
4. 本图检查关系为 `A = 最终白色像素数 × 0.0990358383 × 0.0990358383`。若不满足，先查标定、二值极性和测量选项。

## 13. 在原图上标出主干和棘

1. 在 ShaftPositive Mask 上设阈值 255–255，运行 `Edit > Selection > Create Selection`，按 T 加入 ROI Manager，命名 `Shaft`。
2. 在最终棘 Mask 上做同样操作，命名 `Spines`。它可以是包含多个不相连部分的复合 ROI。
3. 保存全部分析 ROI 为 ZIP，再在用于标注的管理器列表中只保留需要显示的 Shaft、Spines，必要时加上中心线。
4. 用 ROI Manager 的 Properties 给 Shaft 设置绿色描边、Spines 设置品红描边。不要填充遮住原始细节。
5. 激活**相同尺寸、相同原点的原始图**，使用 Show All 或 `Image > Overlay > From ROI Manager`，检查轮廓是否落在对应结构上。

核验：使用裁剪 Mask 得到的 ROI，只能放到相同裁剪原点的图上；不能直接叠到全幅原图。程序的 `source_*` 文件是全幅坐标，带编号的 `01_*`～`11_*` 是裁剪坐标。

## 14. 保存可复核的原图标记

1. 将含 Overlay 的图另存 TIFF，并保存 ROI ZIP，便于以后修改轮廓。
2. 使用 `Image > Overlay > Flatten` 得到一张永久嵌入颜色标记的 RGB 副本，`File > Save As > PNG...` 保存。
3. 若还需要带空间标定的展示图，再保存 RGB TIFF；同时保留原始灰度 TIFF/CZI，后续定量仍使用灰度图与二值 Mask。
4. 至少保存：原图、Centerline/范围/主干/剔除 ROI、最终棘 Mask、面积与长度结果、标记 PNG，以及参数记录。

核验：PNG 上的轮廓在关闭 ROI Manager 后仍可看到；原始 CZI 没有被改写。[Overlay 与 Flatten](https://imagej.net/ij/docs/guide/146-28.html)

## 用菜单核对程序已经导出的数值

1. 打开同一次结果目录的 `source_raw.tif`，在 Properties 核对标定。
2. `File > Import > XY Coordinates...` 选择 `centerline_xy.txt`，再 `Edit > Selection > Area to Line` 转为开放折线；不拟合样条，关闭 Limit to threshold，Measure 得 L。
3. 打开 `source_spine_mask.tif`，Select None，阈值 255–255，勾 Area 和 Limit to threshold、Redirect=None，Measure 得 A。
4. 对照 `results.csv` 的 `spine_area_um2`、`dendrite_length_um` 和 `area_per_length_um2_per_um`。正常只应存在坐标/标定浮点精度差异。
5. 将人工最终 Mask 与程序 Mask 在相同坐标下叠加，检查漏检、误检位置。面积相近不代表像素位置一致，单纯复算程序 Mask 也不能证明生物学分割正确。

[ImageJ 坐标导入与文件保存](https://imagej.net/ij/docs/guide/146-26.html)

## 使用程序涂抹修正后的核验

在程序中先 Analyze，再从 Edit 选择 `Paint erase spines`，在原图上涂掉误识别的分叉；`Paint add spines` 用于补入漏检的阳性区域。`Brush diameter (um)` 控制笔刷直径，松开鼠标后更新面积；`Undo brush` 撤销上一笔，`Clear brush edits` 恢复自动结果。最后运行 Save result 单元格保存。

添加仅纳入测量范围内、主干与排除多边形之外且达到荧光阈值的像素，可人工覆盖自动连接关系、面积和延伸筛选。低于阈值的棘需要先调低荧光阈值。修改阈值会保留涂抹记录；重新选框或清除中心线会清除这些记录。

`source_spine_mask.tif` 是包含涂抹修正的最终掩膜，仍按前面的菜单步骤测量它。另分别测量 `source_auto_spine_mask.tif`、`source_manual_added_mask.tif` 和 `source_manual_removed_mask.tif`，应满足：**最终面积 = 自动面积 + 手动增加面积 − 手动删除面积**。这些修正掩膜表示相对于当前参数自动结果的净变化。`components.csv` 记录自动筛选，`final_components.csv` 对应最终掩膜；`session.json` 保存每一笔在原图上的位置和物理半径。
