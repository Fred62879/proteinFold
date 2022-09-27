from pymol import cmd


pdb_fn='/media/fred/Local Disk/Projects/bioinfo/data/input/idr/pdbs/2AZE.pdb'
pred_fn='/media/fred/Local Disk/Projects/bioinfo/data/output/idr_af_full/poly_g_6_fasta/2AZE.fasta/ranked_0_removed_linker.pdb'

cmd.load(pdb_fn, 'native')
cmd.load(pred_fn, 'pred')

#cmd.select('native_A','native and chain A')
#cmd.select('pred_A','pred and chain A')
#cmd.select('native_B','native and chain B')
#cmd.select('pred_B','pred and chain B')
#cmd.super('native_A','pred_A')
#print(cmd.rms_cur('native_B','pred_B'))

cmd.select('native_AB','native and not chain C')
cmd.select('pred_AB','pred and not chain C')
cmd.select('native_C','native and chain C')
cmd.select('pred_C','pred and chain C')
a = cmd.super('native_AB','pred_AB')
print(a)
print(cmd.rms_cur('native_AB and backbone','pred_AB and backbone
'))
print(cmd.rms_cur('native_C','pred_C'))

cmd.quit()

'''
#arr = os.listdir('../../data/input/dibs/pdbs')
#arr = [ele.split('.')[0].upper() for ele in arr]
#res = reduce(lambda cur, acc: acc + ',' + cur, arr, '')
for name in arr:
    old = join('../../data/input/dibs/pdbs', name.lower()+'.pdb')
    if not exists(old):
        continue
    #new = join('../../data/input/dibs/pdbs', name.split('.')[0][-4:].upper() + '.pdb')
    new = join('../../data/input/dibs/pdbs', name.upper()+'.pdb')
    os.rename(old, new)
'''
