const { allelic_balance_high_quality } = require('./index')
const { high_quality } = require('./index')
const { trio_quality } = require('./index')
const { uniparental_disomy } = require('./index')
const { sample_meets_dominant } = require('./index')
const { dominant } = require('./index')
const { denovo } = require('./index')
const { x_linked_denovo } = require('./index')
const { homozygous_recessive } = require('./index')
const { x_linked_homozygous_recessive } = require('./index')
const { compound_heterozygous_side } = require('./index')
const { remove_compound_heterozygous_side } = require('./index')
const { proband_has_variant } = require('./index')
const { present_in_database } = require('./index')
const { includes_filter } = require('./index')
const { match_filter } = require('./index')
const { gte_filter } = require('./index')
const { lte_filter } = require('./index')
const { maf_filter } = require('./index')

test.each([
  ['allelic_balance_high_quality: homozygous REF with high AB', { alts: 0, AB: 0.2 }, false],
  ['allelic_balance_high_quality: homozygous REF with appropriate AB', { alts: 0, AB: 0.001 }, true],
  ['allelic_balance_high_quality: heterozygous with low AB', { alts: 1, AB: 0.01 }, false],
  ['allelic_balance_high_quality: heterozygous with high AB', { alts: 1, AB: 0.9 }, false],
  ['allelic_balance_high_quality: heterozygous with appropriate AB', { alts: 1, AB: 0.3 }, true],
  ['allelic_balance_high_quality: homozygous REF with appropriate AB', { alts: 2, AB: 0.99 }, true],
])('%s: %j equals %p', (a, b, expected) => {
  expect(allelic_balance_high_quality(b)).toBe(expected)
})

test.each([
  ['high_quality: can set defaults', { DP: 100, GQ: 5 }, 1, 1, true],
  ['high_quality: high DP but low GQ', { DP: 100, GQ: 5 }, undefined, undefined, false],
  ['high_quality: low DP but high GQ', { DP: 1, GQ: 75 }, undefined, undefined, false],
])('%s: %j equals %p. DP: %s GQ: %s', (a, b, depth, gq, expected) => {
  expect(high_quality(b, depth, gq)).toBe(expected)
})

test.each([
  [
    'trio_quality: all have good quality',
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    true,
  ],
  [
    'trio_quality: kid is undefined',
    undefined,
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    false,
  ],
  [
    'trio_quality: kid is {}',
    {},
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    false,
  ],
  [
    'trio_quality: mom is undefined',
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    undefined,
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    true,
  ],
  [
    'trio_quality: dad is undefined',
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    undefined,
    true,
  ],
  [
    'trio_quality: kid has bad quality',
    { DP: 1, GQ: 1, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    false,
  ],
  [
    'trio_quality: mom has bad quality',
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 1, GQ: 1, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    false,
  ],
  [
    'trio_quality: dad has bad quality',
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 100, GQ: 100, alts: 1, AB: 0.5 },
    { DP: 1, GQ: 1, alts: 1, AB: 0.5 },
    false,
  ],
])('%s: kid: %j mom: %j dad: %j equals %p.', (a, kid, mom, dad, expected) => {
  expect(trio_quality(kid, mom, dad)).toBe(expected)
})

test.each([
  ['uniparental_disomy: can handle missing dad', { alts: 2 }, { alts: 2 }, undefined, false],
  ['uniparental_disomy: can handle missing mom', { alts: 2 }, undefined, { alts: 2 }, false],
  ['uniparental_disomy: can handle missing parents', { alts: 2 }, undefined, undefined, false],
  ['uniparental_disomy: Not inherited (homozygous recessive)', { alts: 2 }, { alts: 1 }, { alts: 1 }, false],
  ['uniparental_disomy: 2 ALT alleles inherited from mom', { alts: 2 }, { alts: 2 }, { alts: 0 }, true],
  ['uniparental_disomy: 2 ALT alleles inherited from dad', { alts: 2 }, { alts: 0 }, { alts: 2 }, true],
  ['uniparental_disomy: 2 REF alleles inherited from mom', { alts: 0 }, { alts: 0 }, { alts: 2 }, true],
  ['uniparental_disomy: 2 REF alleles inherited from dad', { alts: 0 }, { alts: 2 }, { alts: 0 }, true],
  ['uniparental_disomy: proband is heterozygous', { alts: 1 }, { alts: 2 }, { alts: 0 }, false],
  [
    'uniparental_disomy: proband is reference and parents are heterozygous',
    { alts: 0 },
    { alts: 1 },
    { alts: 1 },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(uniparental_disomy(kid, mom, dad)).toBe(expected)
})

test.each([
  ['denovo: can handle missing parents', { alts: 1 }, { alts: 0 }, undefined, false],
  ['denovo: can handle missing parents', { alts: 1 }, undefined, undefined, false],
  ['denovo: parents are both homozygous REF', { alts: 1 }, { alts: 0, AD: [1, 0] }, { alts: 0, AD: [1, 0] }, true],
  ['denovo: proband is homozygous', { alts: 2 }, { alts: 0, AD: [1, 0] }, { alts: 0, AD: [1, 0] }, false],
  ['denovo: mom is heterozygous', { alts: 1 }, { alts: 1, AD: [1, 0] }, { alts: 0, AD: [1, 0] }, false],
  ['denovo: dad is heterozygous', { alts: 1 }, { alts: 0, AD: [1, 0] }, { alts: 1, AD: [1, 0] }, false],
  [
    'denovo: parents are both homozygous REF but mom has some AD',
    { alts: 1 },
    { alts: 0, AD: [1, 1] },
    { alts: 0, AD: [1, 0] },
    false,
  ],
  [
    'denovo: parents are both homozygous REF but dad has some AD',
    { alts: 1 },
    { alts: 0, AD: [1, 0] },
    { alts: 0, AD: [1, 1] },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(denovo(kid, mom, dad)).toBe(expected)
})

test.each([
  ['sample_meets_dominant: can handle missing data', {}, false],
  ['sample_meets_dominant: affected and heterozygous', { affected: true, alts: 1 }, true],
  ['sample_meets_dominant: affected and homozygous REF', { affected: true, alts: 0 }, false],
  ['sample_meets_dominant: affected and homozygous ALT', { affected: true, alts: 2 }, false],
  ['sample_meets_dominant: unaffected and homozygous REF', { affected: false, alts: 0 }, true],
  ['sample_meets_dominant: unaffected and heterozygous', { affected: false, alts: 1 }, false],
  ['sample_meets_dominant: unaffected and homozygous ALT', { affected: false, alts: 2 }, false],
])('%s: Data: %j', (a, data, expected) => {
  expect(sample_meets_dominant(data)).toBe(expected)
})

test.each([
  ['dominant: can handle missing dad', { alts: 1 }, { alts: 0, affected: false }, undefined, false],
  ['dominant: can handle missing mom', { alts: 1 }, undefined, { alts: 0, affected: false }, false],
  [
    'dominant: both parents are unaffected',
    { alts: 1 },
    { alts: 0, affected: false },
    { alts: 0, affected: false },
    false,
  ],
  ['dominant: can handle missing both parents', { alts: 1 }, undefined, undefined, false],
  [
    'dominant: proband is homozygous REF',
    { alts: 0 },
    { alts: 0, affected: false },
    { alts: 1, affected: true },
    false,
  ],
  [
    'dominant: proband is homozygous ALT',
    { alts: 2 },
    { alts: 0, affected: false },
    { alts: 1, affected: true },
    false,
  ],
  [
    'dominant: mom is affected and heterozygous and dad is unaffected and REF',
    { alts: 1 },
    { alts: 1, affected: true },
    { alts: 0, affected: false },
    true,
  ],
  [
    'dominant: dad is affected and heterozygous and mom is unaffected and REF',
    { alts: 1 },
    { alts: 0, affected: false },
    { alts: 1, affected: true },
    true,
  ],
  [
    'dominant: mom is unaffected and heterozygous',
    { alts: 1 },
    { alts: 1, affected: false },
    { alts: 0, affected: false },
    false,
  ],
  [
    'dominant: dad is unaffected and heterozygous',
    { alts: 1 },
    { alts: 0, affected: false },
    { alts: 1, affected: false },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(dominant(kid, mom, dad)).toBe(expected)
})

test.each([
  ['x_linked_denovo: can handle missing parents', { alts: 1 }, { alts: 0 }, undefined, true],
  ['x_linked_denovo: can handle missing parents', { alts: 1 }, undefined, undefined, true],
  ['x_linked_denovo: can handle missing sex', { alts: 1 }, { alts: 0 }, { alts: 0 }, true],
  [
    'x_linked_denovo: male proband w/ homozygous(hemi but not in VCF) and parents are both homozygous REF',
    { alts: 2, sex: 'male' },
    { alts: 0 },
    { alts: 0 },
    true,
  ],
  [
    'x_linked_denovo: male proband w/ hemizygous and parents are both homozygous REF',
    { alts: 1, sex: 'male' },
    { alts: 0 },
    { alts: 0 },
    true,
  ],
  [
    'x_linked_denovo: female proband w/ heterozygous and parents are both homozygous REF',
    { alts: 1, sex: 'female' },
    { alts: 0 },
    { alts: 0 },
    true,
  ],
  [
    'x_linked_denovo: female proband w/ homozygous and parents are both homozygous REF',
    { alts: 2, sex: 'female' },
    { alts: 0 },
    { alts: 0 },
    false,
  ],
  ['x_linked_denovo: mom is heterozygous', { alts: 1, sex: 'male' }, { alts: 1 }, { alts: 0 }, false],
  ['x_linked_denovo: dad is hemizygous', { alts: 1, sex: 'male' }, { alts: 0 }, { alts: 1 }, false],
  ['x_linked_denovo: proband is homozygous REF', { alts: 0, sex: 'male' }, { alts: 0 }, { alts: 1 }, false],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(x_linked_denovo(kid, mom, dad)).toBe(expected)
})

test.each([
  ['homozygous_recessive: can handle missing parents', { alts: 2 }, undefined, undefined, true],
  ['homozygous_recessive: mom is heterozygous and dad is heterozygous', { alts: 2 }, { alts: 1 }, { alts: 1 }, true],
  ['homozygous_recessive: mom is homozygous and dad is heterozygous', { alts: 2 }, { alts: 2 }, { alts: 1 }, false],
  ['homozygous_recessive: mom is heterozygous and dad is homozygous', { alts: 2 }, { alts: 1 }, { alts: 2 }, false],
  ['homozygous_recessive: proband is heterozygous', { alts: 1 }, { alts: 1 }, { alts: 2 }, false],
  ['homozygous_recessive: proband is reference', { alts: 0 }, { alts: 1 }, { alts: 2 }, false],
  [
    'homozygous_recessive: proband is homozygous but mom is homozygous REF',
    { alts: 2 },
    { alts: 0 },
    { alts: 2 },
    false,
  ],
  [
    'homozygous_recessive: proband is homozygous but dad is homozygous REF',
    { alts: 2 },
    { alts: 1 },
    { alts: 0 },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(homozygous_recessive(kid, mom, dad)).toBe(expected)
})

test.each([
  [
    'x_linked_homozygous_recessive: can handle missing parents',
    { alts: 1, sex: 'male' },
    { alts: 1, sex: 'female' },
    undefined,
    true,
  ],
  ['x_linked_homozygous_recessive: can handle missing parents', { alts: 1, sex: 'male' }, undefined, undefined, true],
  [
    'x_linked_homozygous_recessive: can handle missing sex',
    { alts: 1 },
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'male' },
    true,
  ],
  [
    'x_linked_homozygous_recessive: hemizygous male and parents are both homozygous REF',
    { alts: 1, sex: 'male' },
    { alts: 0, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: hemizygous (het in VCF) male and parents are both homozygous REF',
    { alts: 2, sex: 'male' },
    { alts: 0, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: heterozygous female and parents are both homozygous REF',
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: homozygous female and parents are both homozygous REF',
    { alts: 2, sex: 'female' },
    { alts: 0, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: hemizygous male and mom is heterozygous',
    { alts: 1, sex: 'male' },
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'male' },
    true,
  ],
  [
    'x_linked_homozygous_recessive: homyzgous REF male and mom is heterozygous',
    { alts: 0, sex: 'male' },
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: hemizygous (het in VCF) male and mom is heterozygous',
    { alts: 2, sex: 'male' },
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'male' },
    true,
  ],
  [
    'x_linked_homozygous_recessive: heterozygous male and mom is heterozygous',
    { alts: 1, sex: 'male' },
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'male' },
    true,
  ],
  [
    'x_linked_homozygous_recessive: homozygous female and both parents are heterozygous',
    { alts: 2, sex: 'female' },
    { alts: 1, sex: 'female' },
    { alts: 1, sex: 'male' },
    true,
  ],
  [
    'x_linked_homozygous_recessive: homozygous REF female and both parents are heterozygous',
    { alts: 0, sex: 'female' },
    { alts: 1, sex: 'female' },
    { alts: 1, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: homozygous female and mom is heterozygous and dad is hemizygous',
    { alts: 2, sex: 'female' },
    { alts: 1, sex: 'female' },
    { alts: 2, sex: 'male' },
    true,
  ],
  [
    'x_linked_homozygous_recessive: heterozygous female and both parents are heterozygous',
    { alts: 1, sex: 'female' },
    { alts: 1, sex: 'female' },
    { alts: 1, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: homozygous female and parents are both homozygous REF',
    { alts: 2, sex: 'female' },
    { alts: 0, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
  [
    'x_linked_homozygous_recessive: homozygous female and mom is heterozygous and dad is homozygous REF',
    { alts: 2, sex: 'female' },
    { alts: 1, sex: 'female' },
    { alts: 0, sex: 'male' },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(x_linked_homozygous_recessive(kid, mom, dad)).toBe(expected)
})

test.each([
  ['compound_heterozygous_side: can handle missing parents', { alts: 1 }, undefined, undefined, true],
  [
    'compound_heterozygous_side: mom is heterozygous and dad is heterozygous',
    { alts: 1 },
    { alts: 1 },
    { alts: 1 },
    true,
  ],
  [
    'compound_heterozygous_side: mom is homozygous and dad is heterozygous',
    { alts: 1 },
    { alts: 2 },
    { alts: 1 },
    false,
  ],
  [
    'compound_heterozygous_side: mom is heterozygous and dad is homozygous',
    { alts: 1 },
    { alts: 1 },
    { alts: 2 },
    false,
  ],
  ['compound_heterozygous_side: mom and dad are both homozygous REF', { alts: 1 }, { alts: 0 }, { alts: 0 }, true],
  ['compound_heterozygous_side: mom is undefined and dad is homozygous REF', { alts: 1 }, undefined, { alts: 0 }, true],
  ['compound_heterozygous_side: dad is undefined and mom is homozygous REF', { alts: 1 }, { alts: 0 }, undefined, true],
  ['compound_heterozygous_side: proband is not heterozygous', { alts: 0 }, { alts: 0 }, undefined, false],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(compound_heterozygous_side(kid, mom, dad)).toBe(expected)
})

test.each([
  ['remove_compound_heterozygous_side: can handle defaults', { compound_heterozygous_side: 'id' }, undefined, false],
  [
    'remove_compound_heterozygous_side: will remove compound_heterozygous_side',
    { compound_heterozygous_side: 'id' },
    undefined,
    false,
  ],
  [
    'remove_compound_heterozygous_side: will remove compound_heterozygous_side',
    { compound_heterozygous_side2: 'id' },
    'compound_heterozygous_side2',
    false,
  ],
  ['remove_compound_heterozygous_side: will not remove compound_heterozygous_side', {}, undefined, true],
])('%s: INFO: %j key: %s equals %p', (a, INFO, key, expected) => {
  expect(remove_compound_heterozygous_side(INFO, key)).toBe(expected)
})

test.each([
  ['proband_has_variant: can handle missing data', {}, false],
  ['proband_has_variant: proband is homozygous REF', { alts: 0 }, false],
  ['proband_has_variant: proband is homozygous ALT', { alts: 2 }, true],
  ['proband_has_variant: proband is heterozygous', { alts: 1 }, true],
])('%s: %j equals %p', (a, info, expected) => {
  expect(proband_has_variant(info)).toBe(expected)
})

test.each([
  ['present_in_database: can handle missing', {}, undefined, undefined, false],
  [
    'present_in_database: can handle default parameters',
    { CSQ: '|x|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    true,
  ],
  ['present_in_database: is not present', { CSQ: '||||||||||||||||||||||||||||||||||||' }, undefined, undefined, false],
  [
    'present_in_database: can handle setting key and position',
    { BSQ: '||x||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    true,
  ],
])('%s: %j key: %s position: %s equals %p', (a, info, key, position, expected) => {
  expect(present_in_database(info, key, position)).toBe(expected)
})

test.each([
  ['includes_filter: can handle missing', {}, undefined, undefined, undefined, false],
  [
    'includes_filter: can handle default parameters',
    { CSQ: '||||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    undefined,
    false,
  ],
  [
    'includes_filter: matches',
    { CSQ: '|missense2|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    ['missense'],
    true,
  ],
  [
    'includes_filter: does not match',
    { CSQ: '|missense2|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    ['frameshift'],
    false,
  ],
  [
    'includes_filter: can handle setting key and position',
    { BSQ: '||frameshift2||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    ['frameshift'],
    true,
  ],
])('%s: %j key: %s position: %s types: %j equals %p', (a, info, key, position, types, expected) => {
  expect(includes_filter(info, key, position, types)).toBe(expected)
})

test.each([
  ['match_filter: can handle missing', {}, undefined, undefined, undefined, false],
  [
    'match_filter: can handle default parameters',
    { CSQ: '||||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    undefined,
    false,
  ],
  [
    'match_filter: matches',
    { CSQ: '|missense|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    ['missense'],
    true,
  ],
  [
    'match_filter: does not match',
    { CSQ: '|frameshift2|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    ['frameshift'],
    false,
  ],
  [
    'match_filter: can handle setting key and position',
    { BSQ: '||frameshift||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    ['frameshift'],
    true,
  ],
])('%s: %j key: %s position: %s types: %j equals %p', (a, info, key, position, types, expected) => {
  expect(match_filter(info, key, position, types)).toBe(expected)
})

test.each([
  ['gte_filter: can handle missing', {}, undefined, undefined, undefined, false],
  [
    'gte_filter: can handle default parameters',
    { CSQ: '||||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    undefined,
    false,
  ],
  ['gte_filter: meets criteria', { CSQ: '|0.1|||||||||||||||||||||||||||||||||||' }, undefined, undefined, 0.01, true],
  [
    'gte_filter: does not meet criteria',
    { CSQ: '|0.01|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    0.1,
    false,
  ],
  [
    'gte_filter: can handle setting key and position',
    { BSQ: '||0.1||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.1,
    true,
  ],
  [
    'gte_filter: can handle position with multiple values',
    { BSQ: '||0.1&0.4||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.3,
    true,
  ],
  [
    'gte_filter: can handle position with multiple values',
    { BSQ: '||0.1&0.4||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.5,
    false,
  ],
])('%s: %j key: %s position: %s cutoff: %s equals %p', (a, info, key, position, cutoff, expected) => {
  expect(gte_filter(info, key, position, cutoff)).toBe(expected)
})

test.each([
  ['lte_filter: can handle missing', {}, undefined, undefined, undefined, false],
  [
    'lte_filter: can handle default parameters',
    { CSQ: '||||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    undefined,
    false,
  ],
  ['lte_filter: meets criteria', { CSQ: '|0.01|||||||||||||||||||||||||||||||||||' }, undefined, undefined, 0.1, true],
  [
    'lte_filter: does not meet criteria',
    { CSQ: '|0.1|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    0.01,
    false,
  ],
  [
    'lte_filter: can handle setting key and position',
    { BSQ: '||0.01||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.1,
    true,
  ],
  [
    'lte_filter: can handle position with multiple values',
    { BSQ: '||0.1&0.4||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.3,
    true,
  ],
  [
    'lte_filter: can handle position with multiple values',
    { BSQ: '||0.6&0.7||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.5,
    false,
  ],
])('%s: %j key: %s position: %s cutoff: %s equals %p', (a, info, key, position, cutoff, expected) => {
  expect(lte_filter(info, key, position, cutoff)).toBe(expected)
})

test.each([
  ['maf_filter: can handle missing', {}, undefined, undefined, undefined, false],
  [
    'maf_filter: can handle default parameters',
    { CSQ: '||||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    undefined,
    true,
  ],
  [
    'maf_filter: can handle missing data',
    { CSQ: '||||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    undefined,
    true,
  ],
  ['maf_filter: meets criteria', { CSQ: '|0.01|||||||||||||||||||||||||||||||||||' }, undefined, undefined, 0.1, true],
  [
    'maf_filter: does not meet criteria',
    { CSQ: '|0.1|||||||||||||||||||||||||||||||||||' },
    undefined,
    undefined,
    0.01,
    false,
  ],
  [
    'maf_filter: can handle setting key and position',
    { BSQ: '||0.01||||||||||||||||||||||||||||||||||' },
    'BSQ',
    2,
    0.1,
    true,
  ],
])('%s: %j key: %s position: %s cutoff: %s equals %p', (a, info, key, position, cutoff, expected) => {
  expect(maf_filter(info, key, position, cutoff)).toBe(expected)
})
