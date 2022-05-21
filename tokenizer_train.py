from tokenizers import Tokenizer
from tokenizers.models import WordLevel

bert_tokenizer = Tokenizer(WordLevel(unk_token="[UNK]"))

from tokenizers import normalizers
from tokenizers.normalizers import NFD

normalizer = normalizers.Sequence([NFD()])
bert_tokenizer.normalizer = normalizer
from tokenizers.pre_tokenizers import Split

bert_tokenizer.pre_tokenizer = Split('\t','removed')

from tokenizers.processors import TemplateProcessing

bert_tokenizer.post_processor = TemplateProcessing(
    single="[CLS] $A [SEP]",
    special_tokens=[
        ("[CLS]", 0),
        ("[SEP]", 1),
	    ("[MSK]", 2),
        ("[PAD]", 3),
        ("[UNK]", 4)
    ],
)
from tokenizers.trainers import WordLevelTrainer

trainer = WordLevelTrainer(
	vocab_size=30000,
	special_tokens=["[CLS]","[SEP]","[UNK]","[MSK]","[PAD]"]
)
files = ['../data/all2.genes']
data = []
for f in files:
	f_h = open(f)
	ls = f_h.readlines()
	for l in ls:
		data.append(l.strip())
bert_tokenizer.train_from_iterator(data, trainer = trainer)

bert_tokenizer.save("../token_model/model_all2/token.json")

from transformers import PreTrainedTokenizerFast
tokenizer = PreTrainedTokenizerFast(tokenizer_file='../token_model/model_all2/token.json')
tokenizer.save_pretrained("../token_model/model_all2/")