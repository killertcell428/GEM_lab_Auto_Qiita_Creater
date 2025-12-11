'use client';

import { useState } from 'react';
import { api, FeedbackCreateRequest } from '@/lib/api';

interface HumanFeedbackPanelProps {
  articleId: string;
  onFeedbackAdded?: () => void;
  defaultPhase?: string;
}

export default function HumanFeedbackPanel({ articleId, onFeedbackAdded, defaultPhase }: HumanFeedbackPanelProps) {
  const [content, setContent] = useState('');
  const [targetSection, setTargetSection] = useState<string>('');
  const [intent, setIntent] = useState<'修正したい' | 'もっと詳しく' | '方針を変えたい' | ''>('');
  const [priority, setPriority] = useState(5);
  const [phase, setPhase] = useState<string>(defaultPhase || '');
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    // バリデーション
    if (!content || !content.trim()) {
      setError('フィードバック内容を入力してください');
      return;
    }
    
    setSubmitting(true);
    setError(null);

    try {
      const request: FeedbackCreateRequest = {
        content: content.trim(),
        target_section: targetSection || undefined,
        intent: intent || undefined,
        priority,
        phase: phase || undefined,
      };
      await api.addFeedback(articleId, request);
      setContent('');
      setTargetSection('');
      setIntent('');
      setPriority(5);
      setPhase(defaultPhase || '');
      if (onFeedbackAdded) {
        onFeedbackAdded();
      }
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'フィードバックの送信に失敗しました';
      setError(errorMessage);
      console.error('フィードバック送信エラー:', err);
    } finally {
      setSubmitting(false);
    }
  };

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
      <h3 className="text-lg font-semibold mb-4">Human Feedbackを追加</h3>
      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      <div className="mb-2">
        <label htmlFor="phase" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
          対象フェーズ
        </label>
        <select
          id="phase"
          value={phase}
          onChange={(e) => setPhase(e.target.value)}
          className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
        >
          <option value="">現在のフェーズ</option>
          <option value="plan">Plan</option>
          <option value="do">Do</option>
          <option value="check">Check</option>
          <option value="act">Act</option>
          <option value="publish">Publish</option>
          <option value="analyze">Analyze</option>
          <option value="other">その他</option>
        </select>
      </div>
      <form onSubmit={handleSubmit} className="space-y-4">
        <div>
          <label htmlFor="feedback-content" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            フィードバック内容 <span className="text-red-500">*</span>
          </label>
          <textarea
            id="feedback-content"
            value={content}
            onChange={(e) => setContent(e.target.value)}
            required
            rows={4}
            className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
            placeholder="例: このセクションをもっと詳しく説明してください"
          />
        </div>

        <div>
          <label htmlFor="target-section" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            対象セクション
          </label>
          <select
            id="target-section"
            value={targetSection}
            onChange={(e) => setTargetSection(e.target.value)}
            className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
          >
            <option value="">AIに任せる</option>
            <option value="全体">全体</option>
            <option value="このセクション">このセクション</option>
          </select>
        </div>

        <div>
          <label htmlFor="intent" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            意図
          </label>
          <select
            id="intent"
            value={intent}
            onChange={(e) => setIntent(e.target.value as any)}
            className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
          >
            <option value="">選択してください</option>
            <option value="修正したい">修正したい</option>
            <option value="もっと詳しく">もっと詳しく</option>
            <option value="方針を変えたい">方針を変えたい</option>
          </select>
        </div>

        <div>
          <label htmlFor="priority" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            優先度: {priority} (1-10、10が最高)
          </label>
          <input
            id="priority"
            type="range"
            min="1"
            max="10"
            value={priority}
            onChange={(e) => setPriority(parseInt(e.target.value))}
            className="w-full"
          />
        </div>

        <button
          type="submit"
          disabled={submitting || !content}
          className="w-full px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
        >
          {submitting ? '送信中...' : 'フィードバックを送信'}
        </button>
      </form>
    </div>
  );
}

