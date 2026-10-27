import pathlib,json,csv,hashlib
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');o=r/'BIT_final_mapping_audit';v=json.loads((o/'BIT_integrity_mapping_final.json').read_text())
row={'resource':'BIT Zenodo14231098 Liver.h5ad / GSE227265','decision':'有条件纳入','url':v['url'],'filename':'Liver.h5ad.complete','bytes':v['bytes'],'sha256':v['sha256'],'completed_download_utc':'2026-10-16T17:36:39.250647+00:00','audit_utc':v['completed_utc'],'cells':v['TF_ID_alignment']['BIT_cells'],'features':v['var']['rows'],'TF_exact_ID_overlap':v['TF_ID_alignment']['intersection'],'TF_same_order':False,'donor_mapping':'未找回','HCC_CCA_mapping':'未找回','original_peak_mapping':'未找回；基因组bin与原峰矩阵不同','archive_fully_verified':False,'allowed_use':'按精确ID关联已有细胞类型与TF活性，作为探索性证据','limits':'不能据barcode后缀推患者；不能按维度套原峰顺序；不能进行患者层面HCC限定推断'}
for suffix,sep in [('csv',','),('tsv','\t')]:
 with (o/f'BIT_resource_audit_final.{suffix}').open('x',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=row,delimiter=sep);w.writeheader();w.writerow(row)
lines=['BIT公开数据目标成员完整性与映射核验报告',v['completed_utc'],'',
'一、文件与下载范围',f"来源：{v['url']}",f"目标：{v['file']}",f"实际大小：{v['bytes']} 字节；与tar成员声明大小一致。",f"SHA256：{v['sha256']}",
'目标成员已完整并可解析；整个2,839,519,229字节归档未完整下载，也未执行整包gzip CRC或归档MD5验证，不称整包完整性通过。下载沿用已验证分块，最后阶段保持非代理。','',
'二、矩阵与注释',f"X维度：{v['matrices']['X']['shape']}；存储条目：{v['matrices']['X']['stored_entries']}。",
f"全量数值与索引扫描：{json.dumps(v['matrices'],ensure_ascii=False)}",f"obs字段：{', '.join(v['obs_fields'])}",f"已有细胞类型计数：{json.dumps(v['celltype_counts'],ensure_ascii=False)}",f"var字段与区间检查：{json.dumps(v['var'],ensure_ascii=False)}",'已有细胞类型属于来源注释，本轮未重新验证生物学正确性；未重聚类或运行研究模型。','',
'三、映射结论',f"与原TF表的精确细胞ID核对：{json.dumps(v['TF_ID_alignment'],ensure_ascii=False)}",'全部18,631个细胞可按精确ID连接，顺序不同，必须使用导出的连接表，不能按位置套用。',
'患者、供体、GSM及HCC/CCA逐细胞映射未找回；不能从barcode后缀推定患者。',
'此对象为6,062,095个基因组区间，原GEO矩阵为1,089,114个峰。未取得两者的明确对应表，因此原峰坐标和顺序缺口仍未解决。var.selected不能自动等同于原峰集合。',
'结论：有条件纳入。可补齐TF表的已有细胞类型标签，用于探索性细胞类型汇总；不支持据此完成患者层面的HCC限定分析或原峰矩阵坐标分析。','',
'四、证据位置',str(o/'BIT_integrity_mapping_final.json'),str(o/'TF_to_BIT_exact_ID_join.tsv'),str(o/'BIT_resource_audit_final.csv'),str(o/'BIT_resource_audit_final.tsv'),
'辅助数据和脚本保留服务器/工作目录。本报告不代表整个研究数据已无缺口；其他队列的既有审计结论保持不变。']
dest=o/'BIT_integrity_mapping_report_20261017.txt'
with dest.open('x',encoding='utf-8') as f:f.write('\n'.join(lines)+'\n')
print(json.dumps({'report':str(dest),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()},ensure_ascii=False))
