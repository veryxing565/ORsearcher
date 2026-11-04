#!/usr/bin/env python3
"""Generate the parameter reference from the installed CLI definitions."""
import argparse
from pathlib import Path

from orsearcher.cli import parser
from orsearcher.presets import PRESETS


def cell(value):
    return str(value).replace('|',r'\|').replace('\n',' ').strip()


def format_value(value):
    if isinstance(value,bool):
        return 'enabled' if value else 'disabled'
    return str(value)


def option_table(actions, command):
    rows = ['| Option | Meaning |', '| --- | --- |']
    for action in actions:
        if not action.option_strings or action.help in (None,argparse.SUPPRESS):
            continue
        flags = ', '.join('`'+flag+'`' for flag in action.option_strings)
        if action.dest in ('replace_hmm','add_hmm'):
            flag = action.option_strings[0]
            flags = f'`{flag} OR FILE`, `{flag} GR FILE`'
        help_text = command._get_formatter()._expand_help(action)
        help_text = (help_text.replace(' (default: None)','').replace(' (default: [])','')
                     .replace('orsearcher_<sample>','`orsearcher_<sample>`'))
        if action.choices:
            help_text += ' Choices: '+', '.join(map(str,action.choices))+'.'
        rows.append(f'| {flags} | {cell(help_text)} |')
    return '\n'.join(rows) + '\n'


def reference():
    commands = parser()._subparsers._group_actions[0].choices
    run = commands['run']
    sections = ['# Command-line reference\n',
        'Run `orsearcher -h` for the current installed options. The tables below list\n'
        'G0 defaults. Lengths are in bases or amino acids as specified. Coverage\n'
        'thresholds use fractions from 0 to 1, and amino-acid identity uses percent.\n']
    for group in run._action_groups:
        if group.title in ('positional arguments','options'):
            continue
        sections += [f'## {group.title}\n',option_table(group._group_actions,run)]
    sections += ['## Preset differences\n',
        'All unlisted settings use G0 values in every preset. An explicitly supplied\n'
        'option overrides the selected preset. Values below refer to one full workflow run.\n']
    by_dest = {action.dest:action for action in run._actions}
    rows = ['| Parameter | G0 | G1 | G2 | G3 |','| --- | --- | --- | --- | --- |']
    for name in PRESETS['G3']:
        action = by_dest[name]
        values = [format_value(preset.get(name,action.default)) for preset in PRESETS.values()]
        rows.append(f'| `{action.option_strings[0]}` | ' + ' | '.join(values) + ' |')
    sections.append('\n'.join(rows)+'\n')
    sections += ['## Other commands\n',
        '`orsearcher check` verifies analysis tools. `orsearcher --version` prints the\n'
        'software version. The remaining commands accept the options below. `refine`\n'
        'and `polish` also use the relevant input and tuning options listed above.\n']
    for name in ('references','prepare-benchmark','benchmark','refine','polish'):
        command = commands[name]
        actions = [a for a in command._actions if a.dest!='help' and (name not in ('refine','polish') or a.dest not in by_dest)]
        sections += [f'### `orsearcher {name}`\n',option_table(actions,command)]
    return '\n'.join(sections)


if __name__=='__main__':
    target = Path(__file__).resolve().parents[1]/'docs/CLI.md'
    target.write_text(reference(),encoding='utf-8',newline='\n')
