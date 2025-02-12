import sys
import yaml
import os
yaml.add_representer(str, lambda dumper, data: dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"'))

tasks_per_node = int(sys.argv[1])
threads_per_task = int(sys.argv[2])


# Execute the code using exec and ensure variables are in the local scope
cpu_bind_masks = list()
bits_per_task = threads_per_task
current_block = 0
num_blocks = 0
current_block_free = 4
for _ in range(tasks_per_node):
    bits_left = bits_per_task
    cpu_bind_masks.append('0x')
    while bits_left >= current_block_free:
        mask = '1' * current_block_free
        cpu_bind_masks[-1] += hex(int(mask, 2))[2:].upper()
        bits_left -= current_block_free
        current_block_free = 4
        num_blocks += 1
    if bits_left > 0:
        mask = '1' * bits_left + '0' * (4 - bits_left)
        cpu_bind_masks[-1] += hex(int(mask, 2))[2:].upper()
        current_block_free -= bits_left
    cpu_bind_masks[-1] += '0' * current_block
    current_block += num_blocks
    num_blocks = 0
result = ",".join(cpu_bind_masks)
print(result)
