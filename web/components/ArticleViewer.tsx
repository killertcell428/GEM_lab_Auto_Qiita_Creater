'use client';

import { useState, useEffect } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { vscDarkPlus } from 'react-syntax-highlighter/dist/esm/styles/prism';

interface ArticleViewerProps {
  markdown: string;
  editable?: boolean;
  onSave?: (content: string) => Promise<void>;
}

export default function ArticleViewer({ markdown, editable = false, onSave }: ArticleViewerProps) {
  const [content, setContent] = useState(markdown);
  const [isEditing, setIsEditing] = useState(false);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [success, setSuccess] = useState(false);

  // markdownが変更されたらcontentを更新
  useEffect(() => {
    if (!isEditing) {
      setContent(markdown);
    }
  }, [markdown, isEditing]);

  const handleSave = async () => {
    if (onSave) {
      setSaving(true);
      setError(null);
      setSuccess(false);
      try {
        await onSave(content);
        setIsEditing(false);
        setSuccess(true);
        setTimeout(() => setSuccess(false), 3000);
      } catch (err) {
        const errorMessage = err instanceof Error ? err.message : '保存に失敗しました';
        setError(errorMessage);
        console.error('保存に失敗しました:', err);
      } finally {
        setSaving(false);
      }
    }
  };

  if (editable && isEditing) {
    return (
      <div>
        {error && (
          <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
            {error}
          </div>
        )}
        {success && (
          <div className="mb-4 p-4 bg-green-100 border border-green-400 text-green-700 rounded">
            保存が完了しました
          </div>
        )}
        <div className="mb-4 flex justify-end space-x-2">
          <button
            onClick={() => {
              setContent(markdown);
              setIsEditing(false);
              setError(null);
              setSuccess(false);
            }}
            disabled={saving}
            className="px-4 py-2 border border-gray-300 rounded-md text-gray-700 dark:text-gray-300 hover:bg-gray-50 dark:hover:bg-gray-700 disabled:opacity-50"
          >
            キャンセル
          </button>
          <button
            onClick={handleSave}
            disabled={saving}
            className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50"
          >
            {saving ? '保存中...' : '保存'}
          </button>
        </div>
        <textarea
          value={content}
          onChange={(e) => {
            setContent(e.target.value);
            setError(null);
          }}
          className="w-full h-96 px-4 py-2 border border-gray-300 rounded-md font-mono text-sm dark:bg-gray-700 dark:border-gray-600 dark:text-white"
          placeholder="Markdown形式で記事を記述してください"
        />
      </div>
    );
  }

  return (
    <div>
      {editable && (
        <div className="mb-4 flex justify-end">
          <button
            onClick={() => setIsEditing(true)}
            className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
          >
            編集
          </button>
        </div>
      )}
      <div className="prose dark:prose-invert max-w-none">
        <ReactMarkdown
          remarkPlugins={[remarkGfm]}
          components={{
            code({ node, inline, className, children, ...props }) {
              const match = /language-(\w+)/.exec(className || '');
              return !inline && match ? (
                <SyntaxHighlighter
                  style={vscDarkPlus}
                  language={match[1]}
                  PreTag="div"
                  {...props}
                >
                  {String(children).replace(/\n$/, '')}
                </SyntaxHighlighter>
              ) : (
                <code className={className} {...props}>
                  {children}
                </code>
              );
            },
          }}
        >
          {markdown}
        </ReactMarkdown>
      </div>
    </div>
  );
}

